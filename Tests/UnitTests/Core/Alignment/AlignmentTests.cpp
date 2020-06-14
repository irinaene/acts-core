// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Alignment/detail/AlignmentEngine.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Fitter/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <cmath>
#include <random>
#include <string>

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {
using SourceLink = MinimalSourceLink;
using Covariance = BoundSymMatrix;

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

///
/// @brief Contruct a telescope-like detector
///
struct TelescopeTrackingGeometry {
  /// Default constructor for the Cubit tracking geometry
  ///
  /// @param gctx the geometry context for this geometry at building time
  TelescopeTrackingGeometry(std::reference_wrapper<const GeometryContext> gctx)
      : geoContext(gctx) {
    using namespace UnitLiterals;

    // Construct the rotation
    double rotationAngle = 90_degree;
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    rotation.col(0) = xPos;
    rotation.col(1) = yPos;
    rotation.col(2) = zPos;

    // Boundaries of the surfaces
    rBounds =
        std::make_shared<const RectangleBounds>(RectangleBounds(0.1_m, 0.1_m));

    // Material of the surfaces
    MaterialProperties matProp(95.7, 465.2, 28.03, 14., 2.32e-3, 0.5_mm);
    surfaceMaterial = std::make_shared<HomogeneousSurfaceMaterial>(matProp);
  }

  ///
  /// Call operator to build the standard cubic tracking geometry
  ///
  std::shared_ptr<const TrackingGeometry> operator()() {
    using namespace UnitLiterals;

    // Set translation vectors
    std::vector<Vector3D> translations;
    translations.reserve(6);
    translations.push_back({-500_mm, 0., 0.});
    translations.push_back({-300_mm, 0., 0.});
    translations.push_back({-100_mm, 0., 0.});
    translations.push_back({100_mm, 0., 0.});
    translations.push_back({300_mm, 0., 0.});
    translations.push_back({500_mm, 0., 0.});

    // Construct layer configs
    std::vector<CuboidVolumeBuilder::LayerConfig> lConfs;
    lConfs.reserve(6);
    unsigned int i;
    for (i = 0; i < translations.size(); i++) {
      CuboidVolumeBuilder::SurfaceConfig sConf;
      sConf.position = translations[i];
      sConf.rotation = rotation;
      sConf.rBounds = rBounds;
      sConf.surMat = surfaceMaterial;
      // The thickness to construct the associated detector element
      sConf.thickness = 1._um;
      sConf.detElementConstructor =
          [](std::shared_ptr<const Transform3D> trans,
             std::shared_ptr<const RectangleBounds> bounds, double thickness) {
            return new Test::DetectorElementStub(trans, bounds, thickness);
          };
      CuboidVolumeBuilder::LayerConfig lConf;
      lConf.surfaceCfg = sConf;
      lConfs.push_back(lConf);
    }

    // Construct volume config
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0., 0., 0.};
    vConf.length = {1.2_m, 1._m, 1._m};
    vConf.layerCfg = lConfs;
    vConf.name = "Tracker";

    // Construct volume builder config
    CuboidVolumeBuilder::Config conf;
    conf.position = {0., 0., 0.};
    conf.length = {1.2_m, 1._m, 1._m};
    conf.volumeCfg = {vConf};  // one volume

    // Build detector
    CuboidVolumeBuilder cvb;
    cvb.setConfig(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto& vb) {
          return cvb.trackingVolume(context, inner, vb);
        });
    TrackingGeometryBuilder tgb(tgbCfg);
    std::shared_ptr<const TrackingGeometry> detector =
        tgb.trackingGeometry(tgContext);

    // Build and return tracking geometry
    return detector;
  }

  RotationMatrix3D rotation = RotationMatrix3D::Identity();
  std::shared_ptr<const RectangleBounds> rBounds = nullptr;
  std::shared_ptr<const ISurfaceMaterial> surfaceMaterial = nullptr;

  std::reference_wrapper<const GeometryContext> geoContext;
};

///
/// @brief Unit test for Kalman fitter with measurements along the x-axis
///
BOOST_AUTO_TEST_CASE(Alignment_zero_field) {
  // Build detector
  TelescopeTrackingGeometry tGeometry(tgContext);
  auto detector = tGeometry();

  // Get the surfaces;
  std::vector<const Surface*> surfaces;
  surfaces.reserve(6);
  detector->visitSurfaces([&](const Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // Create measurements (assuming they are for a linear track parallel to
  // global x-axis)
  std::vector<FittableMeasurement<SourceLink>> measurements;
  measurements.reserve(6);
  Vector2D lPosCenter{10_mm, 10_mm};
  std::array<double, 2> resolution = {30_um, 50_um};
  ActsSymMatrixD<2> cov2D;
  cov2D << resolution[eLOC_0] * resolution[eLOC_0], 0., 0.,
      resolution[eLOC_1] * resolution[eLOC_1];
  for (const auto& surface : surfaces) {
    // 2D measurements
    double dx = resolution[eLOC_0] * gauss(generator);
    double dy = resolution[eLOC_1] * gauss(generator);
    MeasurementType<eLOC_0, eLOC_1> m01(surface->getSharedPtr(), {}, cov2D,
                                        lPosCenter[eLOC_0] + dx,
                                        lPosCenter[eLOC_1] + dy);
    measurements.push_back(std::move(m01));
  }

  // Make a vector of source links as input to the KF
  std::vector<SourceLink> sourcelinks;
  std::transform(measurements.begin(), measurements.end(),
                 std::back_inserter(sourcelinks),
                 [](const auto& m) { return SourceLink{&m}; });

  // The KalmanFitter - we use the eigen stepper for covariance transport
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  // Set initial parameters for the particle track
  Covariance cov;
  cov << std::pow(100_um, 2), 0., 0., 0., 0., 0., 0., std::pow(100_um, 2), 0.,
      0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0.,
      0., 0., 0.01, 0., 0., 0., 0., 0., 0., 1.;

  Vector3D rPos(-1_m, 100_um * gauss(generator), 100_um * gauss(generator));
  Vector3D rMom(1_GeV, 0.025_GeV * gauss(generator),
                0.025_GeV * gauss(generator));

  SingleCurvilinearTrackParameters<ChargedPolicy> rStart(cov, rPos, rMom, 1.,
                                                         42.);

  const Surface* rSurface = &rStart.referenceSurface();

  using Updater = GainMatrixUpdater<BoundParameters>;
  using Smoother = GainMatrixSmoother<BoundParameters>;
  using KalmanFitter = KalmanFitter<RecoPropagator, Updater, Smoother>;

  KalmanFitter kFitter(rPropagator,
                       getDefaultLogger("KalmanFilter", Logging::WARNING));

  KalmanFitterOptions<VoidOutlierFinder> kfOptions(
      tgContext, mfContext, calContext, VoidOutlierFinder(), rSurface);

  // Fit the track
  auto fitRes = kFitter.fit(sourcelinks, rStart, kfOptions);
  BOOST_CHECK(fitRes.ok());
  auto& fittedTrack = *fitRes;

  // Calculate the global track parameters covariance
  const auto& globalTrackParamsCov = detail::globalTrackParametersCovariance(
      fittedTrack.fittedStates, fittedTrack.trackTip);

  // Define the surfaces to be aligned
  std::unordered_map<const Surface*, size_t> alignSurfaces;
  unsigned int iSurface = 0;
  for (const auto& surface : surfaces) {
    // Missing out the forth layer
    if (surface->geoID().layer() == 8) {
      continue;
    }
    alignSurfaces.emplace(surface, iSurface);
    iSurface++;
  }
  // Check the number of aligned surfaces
  BOOST_CHECK_EQUAL(alignSurfaces.size(), 5);

  // Calculate the alignment state for this fitted track
  const auto& alignState = detail::trackAlignmentState(
      fittedTrack.fittedStates, fittedTrack.trackTip, globalTrackParamsCov,
      alignSurfaces);
  std::cout << "Chi2/ndf: " << alignState.chi2 / alignState.alignmentDof
            << std::endl;

  // Check the dimensions
  BOOST_CHECK_EQUAL(alignState.measurementDim, 12);
  BOOST_CHECK_EQUAL(alignState.trackParametersDim, 36);
  // Check the alignment dof
  BOOST_CHECK_EQUAL(alignState.alignmentDof, 30);
  BOOST_CHECK_EQUAL(alignState.alignedSurfaces.size(), 5);
  // Check the measurements covariance
  BOOST_CHECK_EQUAL(alignState.measurementCovariance.rows(), 12);
  const ActsSymMatrixD<2> measCov =
      alignState.measurementCovariance.block<2, 2>(2, 2);
  CHECK_CLOSE_ABS(measCov, cov2D, 1e-10);
  // Check the track parameters covariance matrix. Its rows/columns scales with
  // the number of measurement states
  BOOST_CHECK_EQUAL(alignState.trackParametersCovariance.rows(), 36);
  CHECK_CLOSE_ABS(alignState.trackParametersCovariance,
                  globalTrackParamsCov.first, 1e-10);
  // Check the projection matrix
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.rows(), 12);
  BOOST_CHECK_EQUAL(alignState.projectionMatrix.cols(), 36);
  const ActsMatrixD<2, 6> proj = alignState.projectionMatrix.block<2, 6>(0, 0);
  const ActsMatrixD<2, 6> refProj = ActsMatrixD<2, 6>::Identity();
  CHECK_CLOSE_ABS(proj, refProj, 1e-10);
  // Check the residual
  BOOST_CHECK_EQUAL(alignState.residual.size(), 12);
  // Check the residual covariance
  BOOST_CHECK_EQUAL(alignState.residualCovariance.rows(), 12);
  // Check the chi2 derivative
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2Derivative.size(), 30);
  BOOST_CHECK_EQUAL(alignState.alignmentToChi2SecondDerivative.rows(), 30);
}

}  // namespace Test
}  // namespace Acts
