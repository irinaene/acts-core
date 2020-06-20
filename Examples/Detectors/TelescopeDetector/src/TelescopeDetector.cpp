// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TelescopeDetector/TelescopeDetector.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/Framework/IContextDecorator.hpp"
#include "ACTFW/TelescopeDetector/BuildTelescopeDetector.hpp"
#include "ACTFW/TelescopeDetector/TelescopeDetectorElement.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

void TelescopeDetector::addOptions(
    boost::program_options::options_description& /*opt*/) const {}

auto TelescopeDetector::finalize(
    const boost::program_options::variables_map& /*vm*/,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator)
    -> std::pair<TrackingGeometryPtr, ContextDecorators> {
  // --------------------------------------------------------------------------------
  DetectorElement::ContextType nominalContext;

  /// Return the generic detector
  TrackingGeometryPtr gGeometry = FW::Telescope::buildDetector<DetectorElement>(
      nominalContext, detectorStore, std::move(mdecorator));
  ContextDecorators gContextDeocrators = {};
  // return the pair of geometry and empty decorators
  return std::make_pair<TrackingGeometryPtr, ContextDecorators>(
      std::move(gGeometry), std::move(gContextDeocrators));
}
