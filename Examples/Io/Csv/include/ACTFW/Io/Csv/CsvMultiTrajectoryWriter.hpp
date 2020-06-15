// This file is part of the Acts project.
//
// Copyright (C) 2020 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <fstream>
#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"

namespace FW {
namespace Csv {

/// @class CsvMultiTrajectoryWriter
///
/// Write out the tracks reconstructed using Combinatorial Kalman Filter in
/// comma-separated-value format.
///
/// This writes one file per event into the configured output directory. By
/// default it writes to the current working directory.
/// Files are named using the following schema
///
///     event000000001-CKFtracks.csv
///     event000000002-CKFtracks.csv
///
/// and each line in the file corresponds to one track.
class CsvMultiTrajectoryWriter : public WriterT<TrajectoryContainer> {
 public:
  struct Config {
    std::string inputTrajectories;  ///< Input trajectory collection
    std::string outputDir;          ///< where to place output files
    size_t outputPrecision = 6;     ///< floating point precision
    size_t nMeasurementsMin = 9;    ///< Min number of measurements
  };

  /// constructor
  /// @param cfg is the configuration object
  /// @parm level is the output logging level
  CsvMultiTrajectoryWriter(const Config& cfg,
                           Acts::Logging::Level level = Acts::Logging::INFO);

 protected:
  /// @brief Write method called by the base class
  /// @param [in] context is the algorithm context for consistency
  /// @param [in] tracks is the track collection
  ProcessCode writeT(const FW::AlgorithmContext& context,
                     const TrajectoryContainer& trajectories) final override;

 private:
  Config m_cfg;  //!< Nested configuration struct
};

}  // namespace Csv
}  // namespace FW
