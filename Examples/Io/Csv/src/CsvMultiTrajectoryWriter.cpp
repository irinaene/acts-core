// This file is part of the Acts project.
//
// Copyright (C) 2020 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Csv/CsvMultiTrajectoryWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>
#include <string>

FW::Csv::CsvMultiTrajectoryWriter::CsvMultiTrajectoryWriter(
    const FW::Csv::CsvMultiTrajectoryWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT<TrajectoryContainer>(cfg.inputTrajectories,
                                   "CsvMultiTrajectoryWriter", level),
      m_cfg(cfg) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
  }
}

FW::ProcessCode FW::Csv::CsvMultiTrajectoryWriter::writeT(
    const FW::AlgorithmContext& context,
    const TrajectoryContainer& trajectories) {
  // open per-event file
  std::string path = FW::perEventFilepath(m_cfg.outputDir, "CKFtracks.csv",
                                          context.eventNumber);
  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  // write csv header
  os << "track_id, nMeasurements, nOutliers, nHoles, chi2/ndf, "
        "truthMatchProbability";
  os << '\n';

  // write one line per track
  os << std::setprecision(m_cfg.outputPrecision);
  size_t track_id = 0;
  for (const auto& traj : trajectories) {
    // The trajectory entry indices and the multiTrajectory
    const auto& [trackTips, mj] = traj.trajectory();
    if (trackTips.empty()) {
      ACTS_WARNING("Empty multiTrajectory.");
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
      // Check if the reco track has fitted track parameters
      if (not traj.hasTrackParameters(trackTip)) {
        ACTS_WARNING(
            "No fitted track parameters for trajectory with entry index = "
            << trackTip);
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
        continue;
      }
      // Get the majority particle counts
      size_t nMajorityHits = particleHitCount.front().hitCount;

      // write the track info
      os << track_id << ",";
      os << trajState.nMeasurements << ",";
      os << trajState.nOutliers << ",";
      os << trajState.nHoles << ",";
      os << trajState.chi2Sum / trajState.NDF << ",";
      os << nMajorityHits / trajState.nMeasurements << ",";

      os << '\n';
      track_id++;
    }  // end of one trajectory
  }    // end of multi-trajectory

  return ProcessCode::SUCCESS;
}
