// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ACTFW/EventData/PixelMultiTrajectory.hpp"
#include "ACTFW/EventData/PixelSourceLink.hpp"
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/BareAlgorithm.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

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

#include "Acts/Geometry/TelescopeTrackingGeometry.hpp"

#include <cmath>
#include <random>
#include <string>

using namespace Acts::UnitLiterals;

using SourceLink = FW::PixelSourceLink;
using Covariance = Acts::BoundSymMatrix;

template <ParID_t... params>
using MeasurementType = Acts::Measurement<SourceLink, params...>;

///
/// Function to read in one raw track
///
/// @param inputFile Input file name
/// @param iTrack The track number
std::vector<std::pair<double, double>> readData(const std::string& inputFile,
                                                size_t iTrack) {
  // One element in the vector represent a hit in the track
  std::vector<std::pair<double, double>> protoTrack;

  return protoTrack;
}

///
/// Function to prepare the input tracks for Kalman Fitter
///
std::vector<std::vector<SourceLink>> readInputTrajectories(
    const std::string& inputFile, size_t nTrajectories,
    const std::vector<const Acts::Surface*>& surfaces,
    const std::array<double, 2>& localSigma = {5_um, 5_um}) {
  std::vector<std::vector<SourceLink>> trajectories;
  trajectories.reserve(nTrajectories);
  for (unsigned int iTrack = 0; iTrack < nTrajectories; iTrack++) {
    if (iTrack % 10 == 0) {
      std::cout << "Processing track: " << iTrack << "..." << std::endl;
    }

    // @Todo: call the file reader to read in a map of surface id (0, 1, 2, 3,
    // 4, 5) and 2D local coordinates
    std::vector<std::pair<double, double>> hitLocals =
        readData(inputFile, iTrack);
    assert(hitLocals.size() == 6);

    // setup local covariance
    Acts::ActsSymMatrixD<2> cov2D;
    cov2D << localSigma[0] * localSigma[0], 0., 0.,
        localSigma[1] * localSigma[1];

    // read an input trajectory
    std::vector<SourceLink> sourcelinks;
    sourcelinks.reserve(6);
    for (unsigned int iSurface = 0; iSurface < hitLocals.size(); iSurface++) {
      Acts::Vector2D loc;
      const auto& [locx, locy] = hitLocals.at(iSurface);
      loc << locx, locy;
      // push a hit
      sourcelinks.emplace_back(*surfaces.at(iSurface), 2, loc, cov2D);
    }

    // push the sourcelinks into the trajectory container
    trajectories.push_back(sourcelinks);
  }

  return trajectories;
}

///
///
int main(int argc, char* argv[]) {
  int eventNumber = 1;

  FW::WhiteBoard eventStore(Acts::getDefaultLogger(
      "EventStore#" + std::to_string(eventNumber), Acts::Logging::WARNING));
  // If we ever wanted to run algorithms in parallel, this needs to be
  // changed to Algorithm context copies
  FW::AlgorithmContext algCtx(0, eventNumber, eventStore);

  // Create a test context
  Acts::GeometryContext tgContext = Acts::GeometryContext();
  Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();
  Acts::CalibrationContext calContext = Acts::CalibrationContext();

  std::normal_distribution<double> gauss(0., 1.);
  std::default_random_engine generator(42);

  // Build detector
  Acts::TelescopeTrackingGeometry tGeometry(tgContext);
  auto detector = tGeometry();

  // Get the surfaces;
  std::vector<const Acts::Surface*> surfaces;
  surfaces.reserve(6);
  detector->visitSurfaces([&](const Acts::Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      std::cout << "surface " << surface->geoID() << " placed at: ("
                << surface->center(tgContext).transpose() << " )" << std::endl;
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // Read in the trajectories
  size_t nTrajectories = 1000;
  const auto& inputTrajectories =
      readInputTrajectories("input.txt", nTrajectories, surfaces);

  // The KalmanFitter - we use the eigen stepper for covariance transport
  Acts::Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  Acts::ConstantBField bField(Acts::Vector3D(0., 0., 0.));
  using RecoStepper = Acts::EigenStepper<Acts::ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Acts::Propagator<RecoStepper, Acts::Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  using Updater = Acts::GainMatrixUpdater<Acts::BoundParameters>;
  using Smoother = Acts::GainMatrixSmoother<Acts::BoundParameters>;
  using KalmanFitter = Acts::KalmanFitter<RecoPropagator, Updater, Smoother>;
  using KalmanFitterOptions =
      Acts::KalmanFitterOptions<Acts::VoidOutlierFinder>;

  // Construct the KalmanFitter
  KalmanFitter kFitter(
      rPropagator,
      Acts::getDefaultLogger("KalmanFilter", Acts::Logging::WARNING));

  // Construct the starting track parameter
  // Set initial parameters for the particle track
  Covariance cov;
  cov << std::pow(100_um, 2), 0., 0., 0., 0., 0., 0., std::pow(100_um, 2), 0.,
      0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0.,
      0., 0., 0.01, 0., 0., 0., 0., 0., 0., 1.;

  Acts::Vector3D rPos(-10_um, 10_um * gauss(generator),
                      10_um * gauss(generator));
  Acts::Vector3D rMom(4_GeV, 0.025_GeV * gauss(generator),
                      0.025_GeV * gauss(generator));

  Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(
      cov, rPos, rMom, 1., 42.);
  const Acts::Surface* rSurface = &rStart.referenceSurface();

  // Construct the KalmanFitter options
  KalmanFitterOptions kfOptions(tgContext, mfContext, calContext,
                                Acts::VoidOutlierFinder(), rSurface);

  // Prepare the output data with MultiTrajectory
  std::vector<FW::PixelMultiTrajectory> outputTrajectories;
  outputTrajectories.reserve(nTrajectories);
  for (unsigned int iTrack = 0; iTrack < nTrajectories; iTrack++) {
    const auto& sourcelinks = inputTrajectories.at(iTrack);
    // Fit the track
    auto result = kFitter.fit(sourcelinks, rStart, kfOptions);

    if (result.ok()) {
      // Get the fit output object
      const auto& fitOutput = result.value();
      // The track entry indices container. One element here.
      std::vector<size_t> trackTips;
      trackTips.reserve(1);
      trackTips.emplace_back(fitOutput.trackTip);
      // The fitted parameters container. One element (at most) here.
      FW::IndexedParams indexedParams;
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.value();
        //      ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
        //      ACTS_VERBOSE("  position: " << params.position().transpose());
        //       ACTS_VERBOSE("  momentum: " << params.momentum().transpose());
        // Push the fitted parameters to the container
        indexedParams.emplace(fitOutput.trackTip, std::move(params));
      } else {
        //      ACTS_DEBUG("No fitted paramemeters for track " << itrack);
      }
      // Create a PixelMultiTrajectory
      outputTrajectories.emplace_back(std::move(fitOutput.fittedStates),
                                      std::move(trackTips),
                                      std::move(indexedParams));
    } else {
      //   ACTS_WARNING("Fit failed for track " << itrack << " with error"
      //                                       << result.error());
      // Fit failed, but still create an empty PixelMultiTrajectory
      outputTrajectories.push_back(FW::PixelMultiTrajectory());
    }
  }

  algCtx.eventStore.add(outputTrajectories);

  // write reconstruction performance data
  //  TrackFinderPerformanceWriter::Config perfFinder;
  //  perfFinder.inputParticles = inputParticles;
  //  perfFinder.inputHitParticlesMap = clusterReaderCfg.outputHitParticlesMap;
  //  perfFinder.inputProtoTracks = trackFinderCfg.outputProtoTracks;
  //  perfFinder.outputDir = outputDir;

  return 1;
}
