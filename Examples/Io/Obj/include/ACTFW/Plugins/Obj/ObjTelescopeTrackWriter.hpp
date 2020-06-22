// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace FW {

namespace Obj {

/// @class ObjTelescopeTrackWriter
///
/// Write out the tracks reconstructed using Combinatorial Kalman Filter
/// Writes one file per event with form:
///
///     event000000001-Telescopetracks.obj
///     event000000002-Telescopetracks.obj
///
/// One Thread per write call and hence thread safe
class ObjTelescopeTrackWriter
    : public WriterT<std::vector<PixelMultiTrajectory>> {
 public:
  struct Config {
    std::string inputTrajectories;  ///< input (fitted) trajectories collection
    std::string outputDir;          ///< where to place output files
    double outputScalor = 1.0;      ///< scale output values
    size_t outputPrecision = 6;     ///< floating point precision
    size_t maxNumTracks =
        std::numeric_limits<size_t>::max();  ///< maximum number of tracks to
                                             ///< write
  };

  /// Constructor with arguments
  ///
  /// @param cfg configuration struct
  /// @param level Output logging level
  ObjTelescopeTrackWriter(const Config& cfg,
                          Acts::Logging::Level level = Acts::Logging::INFO);

  /// Virtual destructor
  ~ObjTelescopeTrackWriter() override = default;

  /// End-of-run hook
  ProcessCode endRun() final override;

 private:
  Config m_cfg;  ///!< Internal configuration represenation

 protected:
  /// This implementation holds the actual writing method
  /// and is called by the WriterT<>::write interface
  ProcessCode writeT(
      const AlgorithmContext& context,
      const std::vector<PixelMultiTrajectory>& trackCollection) final override;
};
}  // namespace Obj
}  // namespace FW
