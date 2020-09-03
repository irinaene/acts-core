// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/RiddersPropagator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <limits>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;
using namespace Acts::UnitLiterals;

using MagneticField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<
    MagneticField, Acts::StepperExtensionList<Acts::DenseEnvironmentExtension>>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using RiddersPropagator = Acts::RiddersPropagator<Propagator>;

// absolute parameter tolerances for position, direction, and absolute momentum
constexpr auto epsPos = 10_um;
constexpr auto epsDir = 1_mrad;
constexpr auto epsMom = 5_MeV;
// relative covariance tolerance
constexpr auto epsCov = 0.05;

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;

inline std::shared_ptr<const Acts::TrackingGeometry> makeDetector() {
  using namespace Acts;

  // avoid rebuilding the tracking geometry for every propagator
  static std::shared_ptr<const Acts::TrackingGeometry> detector;
  if (not detector) {
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {0., 0., 0.};
    vConf.length = {4_m, 4_m, 4_m};
    vConf.volumeMaterial = std::make_shared<const HomogeneousVolumeMaterial>(
        Acts::Test::makeBeryllium());
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {0., 0., 0.};
    conf.length = {4_m, 4_m, 4_m};
    CuboidVolumeBuilder cvb(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return cvb.trackingVolume(context, inner, nullptr);
        });
    detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);
  }

  return detector;
}

inline Propagator makePropagator(double bz) {
  using namespace Acts;

  MagneticField magField(Acts::Vector3D(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return Propagator(std::move(stepper), Acts::Navigator(makeDetector()));
}

inline RiddersPropagator makeRiddersPropagator(double bz) {
  using namespace Acts;

  MagneticField magField(Acts::Vector3D(0.0, 0.0, bz));
  Stepper stepper(std::move(magField));
  return RiddersPropagator(std::move(stepper), Acts::Navigator(makeDetector()));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationDenseConstant)

// check that the propagation is reversible and self-consistent

// TODO does not seem to work as-is
BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardBackwardTest<Propagator, Acts::ChargedPolicy,
                         Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, epsPos, epsDir, epsMom);
}

// check that reachable surfaces are correctly reached

// True forward/backward tracks do not work with z cylinders
BOOST_DATA_TEST_CASE(ToCylinderAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest<Propagator, Acts::ChargedPolicy, ZCylinderSurfaceBuilder,
                   Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, ZCylinderSurfaceBuilder(),
      epsPos, epsDir, epsMom);
}

BOOST_DATA_TEST_CASE(ToDisc,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest<Propagator, Acts::ChargedPolicy, DiscSurfaceBuilder,
                   Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, DiscSurfaceBuilder(),
      epsPos, epsDir, epsMom);
}

BOOST_DATA_TEST_CASE(ToPlane,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest<Propagator, Acts::ChargedPolicy, PlaneSurfaceBuilder,
                   Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, PlaneSurfaceBuilder(),
      epsPos, epsDir, epsMom);
}

// True forward/backward tracks do not work with z straws
BOOST_DATA_TEST_CASE(ToStrawAlongZ,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceTest<Propagator, Acts::ChargedPolicy, ZStrawSurfaceBuilder,
                   Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinear(phi, theta, p, q), s, ZStrawSurfaceBuilder(),
      epsPos, epsDir, epsMom);
}

// check covariance transport using the ridders propagator for comparison

BOOST_DATA_TEST_CASE(CovarianceCurvilinear,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runForwardComparisonTest<Propagator, RiddersPropagator, Acts::ChargedPolicy,
                           Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s, epsPos,
      epsDir, epsMom, epsCov);
}

// limit theta to ignore the covariance missmatches at high theta for now
BOOST_DATA_TEST_CASE(CovarianceToCylinderAlongZ,
                     ds::phiWithoutAmbiguity* ds::thetaCentral* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest<Propagator, RiddersPropagator, Acts::ChargedPolicy,
                             ZCylinderSurfaceBuilder,
                             Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZCylinderSurfaceBuilder(), epsPos, epsDir, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(CovarianceToDisc,
                     ds::phi* ds::thetaWithoutBeam* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest<Propagator, RiddersPropagator, Acts::ChargedPolicy,
                             DiscSurfaceBuilder,
                             Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      DiscSurfaceBuilder(), epsPos, epsDir, epsMom, epsCov);
}

BOOST_DATA_TEST_CASE(CovarianceToPlane,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest<Propagator, RiddersPropagator, Acts::ChargedPolicy,
                             PlaneSurfaceBuilder,
                             Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      PlaneSurfaceBuilder(), epsPos, epsDir, epsMom, epsCov);
}

// limit theta to ignore the covariance missmatches at high theta for now
BOOST_DATA_TEST_CASE(CovarianceToStrawAlongZ,
                     ds::phiWithoutAmbiguity* ds::thetaCentral* ds::absMomentum*
                         ds::chargeNonZero* ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, bz) {
  runToSurfaceComparisonTest<Propagator, RiddersPropagator, Acts::ChargedPolicy,
                             ZStrawSurfaceBuilder,
                             Acts::DenseStepperPropagatorOptions>(
      makePropagator(bz), makeRiddersPropagator(bz), geoCtx, magCtx,
      makeParametersCurvilinearWithCovariance(phi, theta, p, q), s,
      ZStrawSurfaceBuilder(), epsPos, epsDir, epsMom, epsCov);
}

BOOST_AUTO_TEST_SUITE_END()