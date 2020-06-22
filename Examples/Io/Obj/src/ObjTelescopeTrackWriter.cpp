// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Obj/ObjTelescopeTrackWriter.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

using Measurement = Acts::Measurement<FW::PixelSourceLink, Acts::ParDef::eLOC_0,
                                      Acts::ParDef::eLOC_1>;

FW::Obj::ObjTelescopeTrackWriter::ObjTelescopeTrackWriter(
    const FW::Obj::ObjTelescopeTrackWriter::Config& cfg,
    Acts::Logging::Level level)
    : WriterT(cfg.inputTrajectories, "ObjTelescopeTrackWriter", level),
      m_cfg(cfg) {
  if (m_cfg.inputTrajectories.empty()) {
    throw std::invalid_argument("Missing input trajectories");
  }
}

FW::ProcessCode FW::Obj::ObjTelescopeTrackWriter::endRun() {
  return FW::ProcessCode::SUCCESS;
}

FW::ProcessCode FW::Obj::ObjTelescopeTrackWriter::writeT(
    const AlgorithmContext& context,
    const std::vector<PixelMultiTrajectory>& trajectories) {
  // open per-event file
  std::string path = FW::perEventFilepath(m_cfg.outputDir, "TelescopeTrack.obj",
                                          context.eventNumber);
  std::ofstream os(path, std::ofstream::out | std::ofstream::trunc);
  if (!os) {
    throw std::ios_base::failure("Could not open '" + path + "' to write");
  }

  ACTS_INFO("Write trajectories to 'TelescopeTrack.obj' in '" << path << "'");

  // Initialize the vertex counter
  unsigned int vCounter = 0;

  size_t iTraj = 0;
  for (const auto& traj : trajectories) {
    // Write only partial tracks if necessary
    if (iTraj == m_cfg.maxNumTracks) {
      break;
    }

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

      // At least three measurements to draw
      if (trajState.nMeasurements > 2) {
        // We start from one
        ++vCounter;

        // visit the states on track
        mj.visitBackwards(trackTip, [&](const auto& state) {
          // we only fill the track states with non-outlier measurement
          auto typeFlags = state.typeFlags();
          if (not typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
            return true;
          }

          // get the measurement
          auto meas = std::get<Measurement>(*state.uncalibrated());

          // get local position
          Acts::Vector2D local(meas.parameters()[Acts::ParDef::eLOC_0],
                               meas.parameters()[Acts::ParDef::eLOC_1]);
          // get global position
          Acts::Vector3D global(0, 0, 0);
          Acts::Vector3D mom(1, 1, 1);
          meas.referenceSurface().localToGlobal(context.geoContext, local, mom,
                                                global);

          // Write the space point
          os << "v " << m_cfg.outputScalor * global.x() << " "
             << m_cfg.outputScalor * global.y() << " "
             << m_cfg.outputScalor * global.z() << '\n';
          return true;
        });
        // Write out the line - only if we have at least two points created
        size_t vBreak = vCounter + trajState.nMeasurements - 1;
        for (; vCounter < vBreak; ++vCounter)
          os << "l " << vCounter << " " << vCounter + 1 << '\n';
      }
    }  // end of one trajectory
    iTraj++;
  }  // end of multi-trajectory

  // return success
  return FW::ProcessCode::SUCCESS;
}
