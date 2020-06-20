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
#include "ACTFW/Io/Performance/TelescopeTrackingPerformanceWriter.hpp"
#include "ACTFW/Io/Root/RootTrajectoryWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/TelescopeDetector/TelescopeDetector.hpp"
#include "ACTFW/TelescopeTracking/TelescopeTrackingAlgorithm.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"

using namespace Acts::UnitLiterals;
using namespace FW;


#include "rapidjson/myrapidjson.h"
///
/// Struct for 2D hit
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
/// Function to read in one raw track
///
std::vector<PixelHit> readTrack(const std::string& fileName, size_t itrack) {
  // One element in the vector represent a 2D hit on the track
  std::vector<PixelHit> track;

  return track;
}

///
/// @brief Struct to read and create the source link tracks
///
struct TelescopeTrackReader {
  /// The detector resolution
  std::array<double, 2> resolution = {5_um, 5_um};

  /// The ordered detector surfaces
  std::vector<const Acts::Surface*> detectorSurfaces;

  /// @param fileName The input file with one line representing one raw track
  /// @param nTracks The number of tracks to process
  ///
  /// @return The created source link tracks
  std::vector<std::vector<PixelSourceLink>> operator()(
      const std::string& fileName, size_t nTracks) const {
    // Create a container for the output tracks
    std::vector<std::vector<PixelSourceLink>> trajectories;
    trajectories.reserve(nTracks);

    // Read in the tracks
    for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
      if (itrack % 100 == 0) {
        std::cout << "Reading in track: " << itrack << "..." << std::endl;
      }

      // @Todo: call the file reader to read in one raw track
      std::vector<PixelHit> hitLocals = readTrack(fileName, itrack);

      // The number of hits should be less or equal to number of provided
      // surfaces?
      assert(hitLocals.size() <= detectorSurfaces.size());

      // Setup local covariance
      Acts::ActsSymMatrixD<2> cov2D;
      cov2D << resolution[0] * resolution[0], 0., 0.,
          resolution[1] * resolution[1];

      // Create the track sourcelinks
      std::vector<PixelSourceLink> sourcelinks;
      sourcelinks.reserve(hitLocals.size());
      for (const auto& hit : hitLocals) {
        Acts::Vector2D loc;
        loc << hit.locX, hit.locY;
        // push a hit
        sourcelinks.emplace_back(*detectorSurfaces.at(hit.surfaceIndex), loc,
                                 cov2D);
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

  // The source link tracks reader
  TelescopeTrackReader trackReader;
  trackReader.detectorSurfaces = surfaces;

  // setup the fitter
  TelescopeTrackingAlgorithm::Config fitter;
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

void json_file_reader {
  
  std::string datafile_name;

  std::FILE* fp = std::fopen(datafile_name.c_str(), "r");
  if(!fp) {
    std::fprintf(stderr, "File opening failed\n");
  }

  char readBuffer[65536];
  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
  rapidjson::Document data;
  data.ParseStream(is);
  if(!data.IsArray()){
    std::frpintf(stderr, "data file is not json array\n");
    exit(2);
  }
  
  rapidjson::Value::ConstValueIterator ev_it = data.Begin();
  rapidjson::Value::ConstValueIterator ev_it_end = data.End();
  while (ev_it != ev_it_end)
  {    
    for (auto& subev : ev_it->GetArray()){
      for (auto& hit : subev["hit_xyz_array"].GetArray()){
        uint16_t xpix   =  hit[0].GetInt();
        uint16_t ypix   =  hit[1].GetInt();
        uint16_t zlayer =  hit[2].GetInt();
      }
    }
    ++ev_it;
  }
}
