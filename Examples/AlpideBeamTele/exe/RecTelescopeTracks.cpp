// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <memory>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Geometry/CommonGeometry.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Utilities/Options.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Utilities/Units.hpp"

#include "TelescopeTrackingPerformanceWriter.hpp"
#include "RootTelescopeTrackWriter.hpp"
#include "ObjTelescopeTrackWriter.hpp"
#include "TelescopeDetector.hpp"
#include "TelescopeTrackingAlgorithm.hpp"

#include "myrapidjson.h"

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

class ClusterPool {
 public:
  void addHit(uint64_t x, uint64_t y, uint64_t z) {
    uint64_t index = x + (y << 16) + (z << 32);
    m_hit_col.push_back(index);
  }

  void buildClusters() {
    std::vector<uint64_t> hit_col_remain = m_hit_col;

    while (!hit_col_remain.empty()) {
      std::vector<uint64_t> hit_col_this_cluster;
      std::vector<uint64_t> hit_col_this_cluster_edge;

      // get first edge seed hit
      // from un-identifed hit to edge hit
      hit_col_this_cluster_edge.push_back(hit_col_remain[0]);
      hit_col_remain.erase(hit_col_remain.begin());

      while (!hit_col_this_cluster_edge.empty()) {
        uint64_t e = hit_col_this_cluster_edge[0];
        uint64_t c = 0x00000001;  // LSB column  x
        uint64_t r = 0x00010000;  // LSB row     y

        //  8 sorround hits search,
        std::vector<uint64_t> sorround_col{e - c + r, e + r,    e + c + r,
                                           e - c,     e + c,    e - c - r,
                                           e - r,     e + c - r};

        for (auto& sr : sorround_col) {
          // only search on un-identifed hits
          auto sr_found_it =
              std::find(hit_col_remain.begin(), hit_col_remain.end(), sr);
          if (sr_found_it != hit_col_remain.end()) {
            // move the found sorround hit
            // from un-identifed hit to an edge hit
            hit_col_this_cluster_edge.push_back(sr);
            hit_col_remain.erase(sr_found_it);
          }
        }

        // after sorround search
        // move from edge hit to cluster hit
        hit_col_this_cluster.push_back(e);
        hit_col_this_cluster_edge.erase(hit_col_this_cluster_edge.begin());
      }

      double cx = 0;
      double cy = 0;
      uint64_t cz = 0;
      for (auto& hit : hit_col_this_cluster) {
        cx += (hit & 0xffff);
        cy += (hit & 0xffff0000) >> 16;
        cz = (hit & 0xffff00000000) >> 32;
      }
      cx /= hit_col_this_cluster.size();
      cy /= hit_col_this_cluster.size();

      m_ccenter_col.push_back(PixelHit{cz, cx, cy});
      m_cluster_col.push_back(std::move(hit_col_this_cluster));
    }
  }

  std::vector<uint64_t> m_hit_col;
  std::vector<std::vector<uint64_t>> m_cluster_col;
  std::vector<PixelHit> m_ccenter_col;
};

///
/// @brief Struct to read and create the source link tracks
///
struct TelescopeTrackReader {
  /// The number of pixels in local x direction
  size_t nPixX = 1024;

  /// The number of pixels in local y direction
  size_t nPixY = 512;

  /// The size of pixel pitch in local x direction
  double pitchX = 26.88_um;

  /// The size of pixel pitch in local y direction
  double pitchY = 29.24_um;

  /// The pixel detector resolution
  std::array<double, 2> resolution = {100_um, 100_um};

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

    std::cout << "There are " << rawTracks.size()
              << " tracks from jsonTrackReader" << std::endl;
    // Loop over the raw track to create the source link track
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
      std::fprintf(stderr, "File %s opening failed\n", fileName.c_str());
      throw std::ios_base::failure("Could not open '" + fileName);
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
      ClusterPool cpool;

      for (auto& subev : ev_it->GetArray()) {
        for (auto& hit : subev["hit_xyz_array"].GetArray()) {
          uint16_t pixX = hit[0].GetInt();
          uint16_t pixY = hit[1].GetInt();
          uint16_t layerIndex = hit[2].GetInt();
          cpool.addHit(pixX, pixY, layerIndex);
          // track.emplace_back(PixelHit{layerIndex,
          //                             (pixX - (nPixX - 1) / 2.) * pitchX,
          //                             (pixY - (nPixY - 1) / 2.) * pitchY});
        }
      }
      cpool.buildClusters();
      ++ev_it;

      std::vector<PixelHit> track = cpool.m_ccenter_col;
      std::cout << "track size " << track.size() << std::endl;

      std::vector<uint16_t> cluster_counter(6, 0);
      for (auto& c : track) {
        cluster_counter[c.surfaceIndex]++;
        c.locX = (c.locX - (nPixX - 1) / 2.) * pitchX;
        c.locY = (c.locY - (nPixY - 1) / 2.) * pitchY;
      }

      bool isgood_data = true;
      for (auto n : cluster_counter) {
        if (n != 1) {
          isgood_data = false;
          break;
        }
      }
      if (!isgood_data) {
        continue;
      }

      std::fprintf(stdout, "\nadd track clusters:");
      for (auto& h : track) {
        std::fprintf(stdout, " [%f, %f, %lu] ", h.locX, h.locY, h.surfaceIndex);
      }
      // push the track to the track container
      rawTracks.push_back(track);
      itrack++;
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
  fitter.inputFileName = inputDir + "/alpide-data.json";
  fitter.outputTrajectories = "trajectories";
  fitter.trackReader = trackReader;
  // The number of tracks you want to process (in default, all of tracks will be
  // read and fitted)
  fitter.maxNumTracks = 10000;
  fitter.fit = TelescopeTrackingAlgorithm::makeFitterFunction(
      trackingGeometry, magneticField, logLevel);
  sequencer.addAlgorithm(
      std::make_shared<TelescopeTrackingAlgorithm>(fitter, logLevel));

  // write tracks as root tree
  RootTelescopeTrackWriter::Config trackRootWriter;
  trackRootWriter.inputTrajectories = fitter.outputTrajectories;
  trackRootWriter.outputDir = outputDir;
  trackRootWriter.outputFilename = "telescope_tracks.root";
  trackRootWriter.outputTreename = "tracks";
  sequencer.addWriter(
      std::make_shared<RootTelescopeTrackWriter>(trackRootWriter, logLevel));

  if (vm["output-obj"].template as<bool>()) {
    // write the tracks (measurements only for the moment) as Csv
    Obj::ObjTelescopeTrackWriter::Config trackObjWriter;
    trackObjWriter.inputTrajectories = fitter.outputTrajectories;
    trackObjWriter.outputDir = outputDir;
    // The number of tracks you want to show (in default, all of tracks will be
    // shown)
    trackObjWriter.maxNumTracks = 100;
    sequencer.addWriter(std::make_shared<Obj::ObjTelescopeTrackWriter>(
        trackObjWriter, logLevel));
  }

  // write reconstruction performance data
  TelescopeTrackingPerformanceWriter::Config perfFitter;
  perfFitter.inputTrajectories = fitter.outputTrajectories;
  perfFitter.outputDir = outputDir;
  sequencer.addWriter(std::make_shared<TelescopeTrackingPerformanceWriter>(
      perfFitter, logLevel));

  return sequencer.run();
}
