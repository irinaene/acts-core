// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Performance/TelescopeTrackingPerformanceWriter.hpp"

#include <TFile.h>
#include <TTree.h>
#include <stdexcept>

#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"
#include "Acts/Utilities/Helpers.hpp"

using Acts::VectorHelpers::eta;

FW::TelescopeTrackingPerformanceWriter::TelescopeTrackingPerformanceWriter(
    FW::TelescopeTrackingPerformanceWriter::Config cfg,
    Acts::Logging::Level lvl)
    : WriterT(cfg.inputTrajectories, "TelescopeTrackingPerformanceWriter", lvl),
      m_cfg(std::move(cfg)),
      m_resPlotTool(m_cfg.resPlotToolConfig, lvl),
      m_effPlotTool(m_cfg.effPlotToolConfig, lvl),
      m_trackSummaryPlotTool(m_cfg.trackSummaryPlotToolConfig, lvl)

{
  // Input track and truth collection name
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories collection");
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

  // initialize the residual and efficiency plots tool
  m_resPlotTool.book(m_resPlotCache);
  m_effPlotTool.book(m_effPlotCache);
  m_trackSummaryPlotTool.book(m_trackSummaryPlotCache);
}

FW::TelescopeTrackingPerformanceWriter::~TelescopeTrackingPerformanceWriter() {
  m_resPlotTool.clear(m_resPlotCache);
  m_effPlotTool.clear(m_effPlotCache);
  m_trackSummaryPlotTool.clear(m_trackSummaryPlotCache);

  if (m_outputFile) {
    m_outputFile->Close();
  }
}

FW::ProcessCode FW::TelescopeTrackingPerformanceWriter::endRun() {
  // fill residual and pull details into additional hists
  m_resPlotTool.refinement(m_resPlotCache);

  if (m_outputFile) {
    m_outputFile->cd();
    m_resPlotTool.write(m_resPlotCache);
    m_effPlotTool.write(m_effPlotCache);
    m_trackSummaryPlotTool.write(m_trackSummaryPlotCache);

    ACTS_INFO("Wrote performance plots to '" << m_outputFile->GetPath() << "'");
  }
  return ProcessCode::SUCCESS;
}

FW::ProcessCode FW::TelescopeTrackingPerformanceWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<PixelMultiTrajectory>& trajectories) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  // Loop over all trajectories
  for (const auto& traj : trajectories) {
    // The trajectory entry indices and the multiTrajectory
    const auto& [trackTips, mj] = traj.trajectory();
    if (trackTips.empty()) {
      ACTS_WARNING("Empty multiTrajectory.");
      continue;
    }

    // Check the size of the trajectory entry indices. For track fitting, there
    // should be at most one trajectory
    if (trackTips.size() > 1) {
      ACTS_ERROR("Track fitting should not result in multiple trajectories.");
      return ProcessCode::ABORT;
    }
    // Get the entry index for the single trajectory
    auto& trackTip = trackTips.front();

    // Select reco track with fitted parameters
    if (not traj.hasTrackParameters(trackTip)) {
      ACTS_WARNING("No fitted track parameters.");
      continue;
    }
    const auto& fittedParameters = traj.trackParameters(trackTip);

    // @Todo: Modify the residual plotting tool and fill the residual plots
    // m_resPlotTool.fill(m_resPlotCache, ctx.geoContext, *ip,
    //                   traj.trackParameters(trackTip));

    // Collect the trajectory summary info
    auto trajState =
        Acts::MultiTrajectoryHelpers::trajectoryState(mj, trackTip);
    // Fill the trajectory summary info
    m_trackSummaryPlotTool.fill(m_trackSummaryPlotCache, fittedParameters,
                                trajState.nStates, trajState.nMeasurements,
                                trajState.nOutliers, trajState.nHoles);
  }

  // Fill the efficiency, how to define

  return ProcessCode::SUCCESS;
}
