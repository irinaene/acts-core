// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <Acts/Utilities/Units.hpp>
#include <memory>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ACTFW/Io/Root/RootTrajectoryWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TelescopeDetector/TelescopeDetector.hpp"
#include "ACTFW/TelescopeTracking/TelescopeTrackingAlgorithm.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::UnitLiterals;
using namespace FW;

///
/// Function to read in one raw track
///
std::vector<std::pair<double, double>> readTrack(const std::string& fileName,
                                                 size_t itrack) {
  // One element in the vector represent a hit in the track
  std::vector<std::pair<double, double>> protoTrack;

  return protoTrack;
}

///
/// @brief Struct to create the source link tracks
///
struct SourceLinkTrackCreator {
  /// The detector resolution
  std::array<double, 2> resolution = {5_um, 5_um};

  /// The ordered detector surfaces
  std::vector<const Acts::Surface*> surfaces;

  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  /// @return The created source link tracks
  std::vector<std::vector<PixelSourceLink>> operator()(
      const std::string& fileName, size_t nTracks) const {
    std::vector<std::vector<PixelSourceLink>> trajectories;
    trajectories.reserve(nTrack);
    for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
      if (itrack % 10 == 0) {
        std::cout << "Processing track: " << itrack << "..." << std::endl;
      }

      // @Todo: call the file reader to read in one raw track
      std::vector<std::pair<double, double>> hitLocals =
          readTrack(fileName, itrack);

      assert(hitLocals.size() == 6);

      // setup local covariance
      Acts::ActsSymMatrixD<2> cov2D;
      cov2D << resolution[0] * resolution[0], 0., 0.,
          resolution[1] * resolution[1];

      // Create the track sourcelinks
      std::vector<PixelSourceLink> sourcelinks;
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
};

int main(int argc, char* argv[]) {
  TelescopeDetector detector;

  // setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc);
  //  detector.addOptions(desc);
  Options::addBFieldOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd =
      std::make_shared<FW::RandomNumbers>(Options::readRandomNumbersConfig(vm));

  // Setup detector geometry
  auto geometry = Geometry::build(vm, detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readBField(vm);

  // Get the surfaces;
  std::vector<const Acts::Surface*> surfaces;
  surfaces.reserve(6);
  trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      surfaces.push_back(surface);
    }
  });
  std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // The source link track creator
  SourceLinkTrackCreator trackCreator;
  trackCreator.surfaces = surfaces;

  // setup the fitter
  TelescopeTrackingAlgorithm::Config fitter;
  fitter.inputFileName = "data.csv";
  fitter.tracksReader = trackCreator;
  fitter.outputTrajectories = "trajectories";
  fitter.randomNumbers = rnd;
  fitter.fit = TelescopeTrackingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<TelescopeTrackingAlgorithm>(fitter, logLevel));

  /*
  // write tracks from fitting
  RootTrajectoryWriter::Config trackWriter;
  trackWriter.inputParticles = inputParticles;
  trackWriter.inputTrajectories = fitter.outputTrajectories;
  trackWriter.outputDir = outputDir;
  trackWriter.outputFilename = "tracks.root";
  trackWriter.outputTreename = "tracks";
  sequencer.addWriter(
      std::make_shared<RootTrajectoryWriter>(trackWriter, logLevel));

  // write reconstruction performance data
  TrackFitterPerformanceWriter::Config perfFitter;
  perfFitter.inputParticles = inputParticles;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.outputDir = outputDir;
  sequencer.addWriter(
      std::make_shared<TrackFitterPerformanceWriter>(perfFitter, logLevel));
*/

  return sequencer.run();
}
