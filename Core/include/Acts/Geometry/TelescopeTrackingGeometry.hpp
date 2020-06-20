// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Alignment/Alignment.hpp"
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
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "Acts/Geometry/DetectorElementStub.hpp"

#include <cmath>
#include <random>
#include <string>

using namespace Acts::UnitLiterals;

namespace Acts {

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

    // Boundaries of the surfaces (ALPIDE SIZE: 13.76256 mm * 30.81896 mm)
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
            return new DetectorElementStub(trans, bounds, thickness);
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

}  // namespace Acts
