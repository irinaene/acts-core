// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include <Acts/Utilities/Units.hpp>
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Io/Performance/TelescopeTrackingPerformanceWriter.hpp"
#include "ACTFW/Io/Root/RootTrajectoryWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TelescopeDetector/TelescopeDetector.hpp"
#include "ACTFW/TelescopeTracking/TelescopeTrackingAlgorithm.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"

#include "rapidjson/myrapidjson.h"

using namespace Acts::UnitLiterals;
using namespace FW;

///
/// Struct for 2D pixel hit
///
struct PixelHit {
  // The correspoinding hit index
  size_t surfaceIndex = 0;

  // The local x coordinate
  double locX = 0;

  // The local y cooridnate
  double locY = 0;
};

///
/// @brief Struct to read and create the source link tracks
///
struct TelescopeTrackReader {
  /// The number of pixels in local x direction
  size_t nPixX = 1054;

  /// The number of pixels in local y direction
  size_t nPixY = 512;

  /// The size of pixel pitch in local x direction
  double pitchX = 26.88_um;

  /// The size of pixel pitch in local y direction
  double pitchY = 29.24_um;

  /// The pixel detector resolution
  std::array<double, 2> resolution = {5_um, 5_um};

  /// The ordered detector surfaces
  std::vector<const Acts::Surface*> detectorSurfaces;

  /// Function to read and create source link tracks as input of fitter
  ///
  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  ///
  /// @return The created source link tracks
  std::vector<std::vector<PixelSourceLink>> operator()(
      const std::string& fileName, size_t nTracks) const {
    // Create a container for the output tracks
    std::vector<std::vector<PixelSourceLink>> sourcelinkTracks;
    sourcelinkTracks.reserve(nTracks);

    // Read in the raw tracks
    std::vector<std::vector<PixelHit>> rawTracks =
        jsonTrackReader(fileName, nTracks);

    // Create the source link tracks with raw tracks
    for (const auto& rtrack : rawTracks) {
      // The number of hits should be less or equal to number of provided
      // surfaces?
      assert(rtrack.size() <= detectorSurfaces.size());

      // Setup local covariance
      Acts::ActsSymMatrixD<2> cov2D;
      cov2D << resolution[0] * resolution[0], 0., 0.,
          resolution[1] * resolution[1];

      // Create the track sourcelinks
      std::vector<PixelSourceLink> sourcelinks;
      sourcelinks.reserve(rtrack.size());
      for (const auto& hit : rtrack) {
        Acts::Vector2D loc;
        loc << hit.locX, hit.locY;
        // push a hit
        sourcelinks.emplace_back(*detectorSurfaces.at(hit.surfaceIndex), loc,
                                 cov2D);
      }
      // push the sourcelinks into the trajectory container
      sourcelinkTracks.push_back(sourcelinks);
    }
    return sourcelinkTracks;
  }

 private:
  /// Function to read and create raw tracks from a json file
  ///
  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  ///
  /// @return The created raw tracks
  std::vector<std::vector<PixelHit>> jsonTrackReader(
      const std::string& fileName, size_t nTracks) const {
    std::FILE* fp = std::fopen(fileName.c_str(), "r");
    if (!fp) {
      std::fprintf(stderr, "File opening failed\n");
    }

    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document data;
    data.ParseStream(is);
    if (!data.IsArray()) {
      std::fprintf(stderr, "data file is not json array\n");
      exit(2);
    }

    rapidjson::Value::ConstValueIterator ev_it = data.Begin();
    rapidjson::Value::ConstValueIterator ev_it_end = data.End();

    std::vector<std::vector<PixelHit>> rawTracks;
    rawTracks.reserve(nTracks);
    size_t itrack = 0;
    while (ev_it != ev_it_end and itrack < nTracks) {
      std::vector<PixelHit> track;
      for (auto& subev : ev_it->GetArray()) {
        for (auto& hit : subev["hit_xyz_array"].GetArray()) {
          uint16_t pixX = hit[0].GetInt();
          uint16_t pixY = hit[1].GetInt();
          uint16_t layerIndex = hit[2].GetInt();
          track.emplace_back(PixelHit{layerIndex, (pixX - nPixX / 2.) * pitchX,
                                      (pixY - nPixY / 2.) * pitchY});
        }
      }
      // push the track to the track container
      rawTracks.push_back(track);
      itrack++;
      ++ev_it;
    }
    return rawTracks;
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

  // The source link tracks reader
  TelescopeTrackReader trackReader;
  trackReader.detectorSurfaces = surfaces;

  // setup the fitter
  TelescopeTrackingAlgorithm::Config fitter;
  //@Todo: add run number information in the file name
  fitter.inputFileName = inputDir + "/data.json";
  fitter.trackReader = trackReader;
  fitter.outputTrajectories = "trajectories";
  fitter.randomNumbers = rnd;
  fitter.fit = TelescopeTrackingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<TelescopeTrackingAlgorithm>(fitter, logLevel));

  // write reconstruction performance data
  TelescopeTrackingPerformanceWriter::Config perfFitter;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.outputDir = outputDir;
  sequencer.addWriter(std::make_shared<TelescopeTrackingPerformanceWriter>(
      perfFitter, logLevel));

  return sequencer.run();
}
