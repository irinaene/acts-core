// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Io/Performance/CKFPerformanceWriterOptions.hpp"

#include <string>
#include <boost/program_options.hpp>

void FW::Options::addCKFPerformanceWriterOptions(
    FW::Options::Description& desc) {
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("use-ml-track-classification", value<bool>()->default_value(false),
      "Use trained neural network to classify reconstructed tracks as "
      "good/duplicate/fake");
}

FW::CKFPerformanceWriter::Config FW::Options::readCKFPerformanceWriterConfig(
    const FW::Options::Variables& variables) {
  using Config = typename FW::CKFPerformanceWriter::Config;
  Config ckfPerfWriterCfg;

  ckfPerfWriterCfg.useMLTrackClassifier =
      variables["use-ml-track-classification"].template as<bool>();

  return ckfPerfWriterCfg;
}