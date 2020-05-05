// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <numeric>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Io/Performance/CKFPerformanceWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

FW::CKFPerformanceWriter::CKFPerformanceWriter(
    FW::CKFPerformanceWriter::Config cfg, Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "CKFPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_fakeRatePlotTool(m_cfg.fakeRatePlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl) {
  // Input track and truth collection name
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
  if (m_cfg.outputFilename.empty()) {
    throw std::invalid_argument("Missing output filename");
  }

  // the output file can not be given externally since TFile accesses to the
  // same file from multiple threads are unsafe.
  // must always be opened internally
  auto path = joinPaths(m_cfg.outputDir, m_cfg.outputFilename);
  m_outputFile = TFile::Open(path.c_str(), "RECREATE");
  if (not m_outputFile) {
    throw std::invalid_argument("Could not open '" + path + "'");
  }

  // initialize the plot tools
  m_effPlotTool.book(m_effPlotCache);
  m_fakeRatePlotTool.book(m_fakeRatePlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

FW::CKFPerformanceWriter::~CKFPerformanceWriter() {
  m_effPlotTool.clear(m_effPlotCache);
  m_fakeRatePlotTool.clear(m_fakeRatePlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);
  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::CKFPerformanceWriter::endRun() {
  if (m_outputFile) {
    m_outputFile->cd();
    m_effPlotTool.write(m_effPlotCache);
    m_fakeRatePlotTool.write(m_fakeRatePlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);
    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::CKFPerformanceWriter::writeT(
    const AlgorithmContext& ctx, const TrajectoryContainer& trajectories) {
  // Read truth particles from input collection
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Counter of truth-matched reco tracks
  std::map<ActsFatras::Barcode, size_t> matched;
  std::map<ActsFatras::Barcode, size_t> unmatched;

  // Loop over all trajectories
  for (const auto& traj : trajectories) {
    if (not traj.hasTrajectory()) {
      ACTS_WARNING("No multiTrajectory available.");
      continue;
    }

    // The trajectory entry indices and the multiTrajectory
    const auto& [trackTips, mj] = traj.trajectory();
    if (trackTips.empty()) {
      ACTS_WARNING("No trajectory entry index found.");
      continue;
    }

    // Loop over all trajectories in a multiTrajectory
    for (const size_t& trackTip : trackTips) {
      // Collect the trajectory summary info
      auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);

      // Reco track selection
      //@TODO: add interface for applying others cuts on reco tracks:
      // -> pT, d0, z0, detector-specific hits/holes number cut
      if (trajState.nMeasurements < m_cfg.nMeasurementsMin) {
        continue;
      }

      // Get the majority truth particle to this track
      std::vector<ParticleHitCount> particleHitCount =
          traj.identifyMajorityParticle(trackTip);
      if (particleHitCount.empty()) {
        ACTS_WARNING(
            "No truth particle associated with this trajectory with entry "
            "index = "
            << trackTip);
      }
      // Find the truth particle with the majority barcode
      // @TODO: pass reconstructed track parameters instead of truth particle to
      // plotting tools.
      ActsFatras::Barcode majorityParticleId =
          particleHitCount.front().particleId;
      auto ip = particles.find(majorityParticleId);
      if (ip == particles.end()) {
        ACTS_WARNING("Majority particle with particleId = "
                     << majorityParticleId
                     << " not found in the particles collection.");
        // Currently, we could only skip those trajectory with
        // majority particle not found (it happens since there is a
        // pre-selection of the truth particles). But this should be fixed.
        continue;
      }

      // Fill the trajectory summary info
      m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, *ip,
                                  trajState.nStates, trajState.nMeasurements,
                                  trajState.nOutliers, trajState.nHoles);

      // Check if the trajectory is matched with truth
      bool isFake = false;
      size_t nMajorityHits = particleHitCount.front().hitCount;
      // Selection of the tracks
      if (nMajorityHits * 1. / trajState.nMeasurements >=
          m_cfg.truthMatchProbMin) {
        matched[majorityParticleId]++;
      } else {
        isFake = true;
        unmatched[majorityParticleId]++;
      }
      // Fill fake rate plots
      m_fakeRatePlotTool.fill(m_fakeRatePlotCache, *ip, isFake);
    }  // end all trajectories in a multiTrajectory
  }    // end all multiTrajectories

  // Loop over all truth particles for efficiency plots and reco details
  // @TODO: add duplication plots
  for (const auto& particle : particles) {
    auto particleId = particle.particleId();
    // Investigate the truth-matched tracks
    size_t nMatchedTracks = 0;
    bool isReconstructed = false;
    auto imatched = matched.find(particleId);
    if (imatched != matched.end()) {
      nMatchedTracks = imatched->second;
      isReconstructed = true;
    }
    // Fill efficiency plots
    m_effPlotTool.fill(m_effPlotCache, particle, isReconstructed);

    // Investigate the fake (i.e. truth-unmatched) tracks
    size_t nFakeTracks = 0;
    auto ifake = unmatched.find(particleId);
    if (ifake != unmatched.end()) {
      nFakeTracks = ifake->second;
    }
    // Fill number of reconstructed/truth-matched/fake tracks for this particle
    m_fakeRatePlotTool.fill(m_fakeRatePlotCache, particle, nMatchedTracks,
                            nFakeTracks);

  }  // end all truth particles

  return ProcessCode::SUCCESS;
}