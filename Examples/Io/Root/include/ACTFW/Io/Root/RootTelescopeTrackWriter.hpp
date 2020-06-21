// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <mutex>
#include <vector>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

class TFile;
class TTree;

namespace FW {

/// @class RootTelescopeTrackWriter
///
/// Write out a trajectory (i.e. a vector of
/// trackState at the moment) into a TTree
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
///
/// Each entry in the TTree corresponds to one trajectory for optimum
/// writing speed. The event number is part of the written data.
///
/// A common file can be provided for to the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing
/// file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class RootTelescopeTrackWriter final
    : public WriterT<std::vector<PixelMultiTrajectory>> {
 public:
  /// @brief The nested configuration struct
  struct Config {
    std::string inputTrajectories;  ///< input (fitted) trajectories collection
    std::string outputDir;          ///< output directory
    std::string outputFilename = "telescope_tracks.root";  ///< output filename
    std::string outputTreename = "tracks";  ///< name of the output tree
    std::string fileMode = "RECREATE";      ///< file access mode
    TFile* rootFile = nullptr;              ///< common root file
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param level Message level declaration
  RootTelescopeTrackWriter(const Config& cfg, Acts::Logging::Level lvl);

  /// Virtual destructor
  ~RootTelescopeTrackWriter() final override;

  /// End-of-run hook
  ProcessCode endRun() final override;

 protected:
  /// @brief Write method called by the base class
  /// @param [in] ctx is the algorithm context for event information
  /// @param [in] trajectories are what to be written out
  ProcessCode writeT(
      const AlgorithmContext& ctx,
      const std::vector<PixelMultiTrajectory>& trajectories) final override;

 private:
  Config m_cfg;             ///< The config class
  std::mutex m_writeMutex;  ///< Mutex used to protect multi-threaded writes
  TFile* m_outputFile{nullptr};  ///< The output file
  TTree* m_outputTree{nullptr};  ///< The output tree
  int m_eventNr{0};              ///< the event number
  int m_trajNr{0};               ///< the trajectory number

  int m_nStates{0};                 ///< number of all states
  int m_nMeasurements{0};           ///< number of states with measurements
  std::vector<int> m_volumeID;      ///< volume identifier
  std::vector<int> m_layerID;       ///< layer identifier
  std::vector<int> m_moduleID;      ///< surface identifier
  std::vector<float> m_lx_hit;      ///< uncalibrated measurement local x
  std::vector<float> m_ly_hit;      ///< uncalibrated measurement local y
  std::vector<float> m_x_hit;       ///< uncalibrated measurement global x
  std::vector<float> m_y_hit;       ///< uncalibrated measurement global y
  std::vector<float> m_z_hit;       ///< uncalibrated measurement global z
  std::vector<float> m_res_x_hit;   ///< hit residual x
  std::vector<float> m_res_y_hit;   ///< hit residual y
  std::vector<float> m_err_x_hit;   ///< hit err x
  std::vector<float> m_err_y_hit;   ///< hit err y
  std::vector<float> m_pull_x_hit;  ///< hit pull x
  std::vector<float> m_pull_y_hit;  ///< hit pull y
  std::vector<int> m_dim_hit;       ///< dimension of measurement

  int m_nPredicted{0};      ///< number of states with predicted parameter
  std::vector<bool> m_prt;  ///< predicted status
  std::vector<float> m_eLOC0_prt;       ///< predicted parameter eLOC0
  std::vector<float> m_eLOC1_prt;       ///< predicted parameter eLOC1
  std::vector<float> m_ePHI_prt;        ///< predicted parameter ePHI
  std::vector<float> m_eTHETA_prt;      ///< predicted parameter eTHETA
  std::vector<float> m_eQOP_prt;        ///< predicted parameter eQOP
  std::vector<float> m_eT_prt;          ///< predicted parameter eT
  std::vector<float> m_res_eLOC0_prt;   ///< predicted parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_prt;   ///< predicted parameter eLOC1 residual
  std::vector<float> m_err_eLOC0_prt;   ///< predicted parameter eLOC0 error
  std::vector<float> m_err_eLOC1_prt;   ///< predicted parameter eLOC1 error
  std::vector<float> m_err_ePHI_prt;    ///< predicted parameter ePHI error
  std::vector<float> m_err_eTHETA_prt;  ///< predicted parameter eTHETA error
  std::vector<float> m_err_eQOP_prt;    ///< predicted parameter eQOP error
  std::vector<float> m_err_eT_prt;      ///< predicted parameter eT error
  std::vector<float> m_pull_eLOC0_prt;  ///< predicted parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_prt;  ///< predicted parameter eLOC1 pull
  std::vector<float> m_x_prt;           ///< predicted global x
  std::vector<float> m_y_prt;           ///< predicted global y
  std::vector<float> m_z_prt;           ///< predicted global z
  std::vector<float> m_px_prt;          ///< predicted momentum px
  std::vector<float> m_py_prt;          ///< predicted momentum py
  std::vector<float> m_pz_prt;          ///< predicted momentum pz
  std::vector<float> m_eta_prt;         ///< predicted momentum eta
  std::vector<float> m_pT_prt;          ///< predicted momentum pT

  int m_nFiltered{0};              ///< number of states with filtered parameter
  std::vector<bool> m_flt;         ///< filtered status
  std::vector<float> m_eLOC0_flt;  ///< filtered parameter eLOC0
  std::vector<float> m_eLOC1_flt;  ///< filtered parameter eLOC1
  std::vector<float> m_ePHI_flt;   ///< filtered parameter ePHI
  std::vector<float> m_eTHETA_flt;      ///< filtered parameter eTHETA
  std::vector<float> m_eQOP_flt;        ///< filtered parameter eQOP
  std::vector<float> m_eT_flt;          ///< filtered parameter eT
  std::vector<float> m_res_eLOC0_flt;   ///< filtered parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_flt;   ///< filtered parameter eLOC1 residual
  std::vector<float> m_err_eLOC0_flt;   ///< filtered parameter eLOC0 error
  std::vector<float> m_err_eLOC1_flt;   ///< filtered parameter eLOC1 error
  std::vector<float> m_err_ePHI_flt;    ///< filtered parameter ePHI error
  std::vector<float> m_err_eTHETA_flt;  ///< filtered parameter eTHETA error
  std::vector<float> m_err_eQOP_flt;    ///< filtered parameter eQOP error
  std::vector<float> m_err_eT_flt;      ///< filtered parameter eT error
  std::vector<float> m_pull_eLOC0_flt;  ///< filtered parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_flt;  ///< filtered parameter eQOP pull
  std::vector<float> m_x_flt;           ///< filtered global x
  std::vector<float> m_y_flt;           ///< filtered global y
  std::vector<float> m_z_flt;           ///< filtered global z
  std::vector<float> m_px_flt;          ///< filtered momentum px
  std::vector<float> m_py_flt;          ///< filtered momentum py
  std::vector<float> m_pz_flt;          ///< filtered momentum pz
  std::vector<float> m_eta_flt;         ///< filtered momentum eta
  std::vector<float> m_pT_flt;          ///< filtered momentum pT
  std::vector<float> m_chi2;            ///< chisq from filtering

  int m_nSmoothed{0};              ///< number of states with smoothed parameter
  std::vector<bool> m_smt;         ///< smoothed status
  std::vector<float> m_eLOC0_smt;  ///< smoothed parameter eLOC0
  std::vector<float> m_eLOC1_smt;  ///< smoothed parameter eLOC1
  std::vector<float> m_ePHI_smt;   ///< smoothed parameter ePHI
  std::vector<float> m_eTHETA_smt;      ///< smoothed parameter eTHETA
  std::vector<float> m_eQOP_smt;        ///< smoothed parameter eQOP
  std::vector<float> m_eT_smt;          ///< smoothed parameter eT
  std::vector<float> m_res_eLOC0_smt;   ///< smoothed parameter eLOC0 residual
  std::vector<float> m_res_eLOC1_smt;   ///< smoothed parameter eLOC1 residual
  std::vector<float> m_err_eLOC0_smt;   ///< smoothed parameter eLOC0 error
  std::vector<float> m_err_eLOC1_smt;   ///< smoothed parameter eLOC1 error
  std::vector<float> m_err_ePHI_smt;    ///< smoothed parameter ePHI error
  std::vector<float> m_err_eTHETA_smt;  ///< smoothed parameter eTHETA error
  std::vector<float> m_err_eQOP_smt;    ///< smoothed parameter eQOP error
  std::vector<float> m_err_eT_smt;      ///< smoothed parameter eT error
  std::vector<float> m_pull_eLOC0_smt;  ///< smoothed parameter eLOC0 pull
  std::vector<float> m_pull_eLOC1_smt;  ///< smoothed parameter eLOC1 pull
  std::vector<float> m_x_smt;           ///< smoothed global x
  std::vector<float> m_y_smt;           ///< smoothed global y
  std::vector<float> m_z_smt;           ///< smoothed global z
  std::vector<float> m_px_smt;          ///< smoothed momentum px
  std::vector<float> m_py_smt;          ///< smoothed momentum py
  std::vector<float> m_pz_smt;          ///< smoothed momentum pz
  std::vector<float> m_eta_smt;         ///< smoothed momentum eta
  std::vector<float> m_pT_smt;          ///< smoothed momentum pT
};

}  // namespace FW
