// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/TelescopeTracking/TelescopeTrackingAlgorithm.hpp"

#include <stdexcept>

#include "ACTFW/EventData/Track.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

FW::TelescopeTrackingAlgorithm::TelescopeTrackingAlgorithm(
    Config cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TelescopeTrackingAlgorithm", level),
      m_cfg(std::move(cfg)) {
  if (m_cfg.inputFileName.empty()) {
    throw std::invalid_argument("Missing input data collection");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

FW::ProcessCode FW::TelescopeTrackingAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(ctx);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  // Read input data
  const std::vector<SourceLinkTrack> sourcelinkTracks =
      m_cfg.trackReader(m_cfg.inputFileName, m_cfg.nTracks);

  // Prepare the output data with MultiTrajectory
  std::vector<PixelMultiTrajectory> trajectories;
  trajectories.reserve(sourcelinkTracks.size());

  for (std::size_t itrack = 0; itrack < sourcelinkTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters
    const auto& trackSourcelinks = sourcelinkTracks[itrack];

    // We can have empty tracks which must give empty fit results
    if (trackSourcelinks.empty()) {
      trajectories.push_back(PixelMultiTrajectory());
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    // Set initial parameters for the particle track
    Acts::BoundSymMatrix cov;
    cov << std::pow(100_um, 2), 0., 0., 0., 0., 0., 0., std::pow(100_um, 2), 0.,
        0., 0., 0., 0., 0., 0.025, 0., 0., 0., 0., 0., 0., 0.025, 0., 0., 0.,
        0., 0., 0., 0.01, 0., 0., 0., 0., 0., 0., 1.;

    Acts::Vector3D rPos(-10_um, 10_um * stdNormal(rng), 10_um * stdNormal(rng));
    Acts::Vector3D rMom(4_GeV, 0.025_GeV * stdNormal(rng),
                        0.025_GeV * stdNormal(rng));

    Acts::SingleCurvilinearTrackParameters<Acts::ChargedPolicy> rStart(
        cov, rPos, rMom, 1., 42.);
    const Acts::Surface* rSurface = &rStart.referenceSurface();

    // Set the KalmanFitter options
    Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions(
        ctx.geoContext, ctx.magFieldContext, ctx.calibContext,
        Acts::VoidOutlierFinder(), rSurface);

    ACTS_DEBUG("Invoke fitter");
    auto result = m_cfg.fit(trackSourcelinks, rStart, kfOptions);
    if (result.ok()) {
      // Get the fit output object
      const auto& fitOutput = result.value();
      // The track entry indices container. One element here.
      std::vector<size_t> trackTips;
      trackTips.reserve(1);
      trackTips.emplace_back(fitOutput.trackTip);
      // The fitted parameters container. One element (at most) here.
      IndexedParams indexedParams;
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.value();
        ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
        ACTS_VERBOSE("  position: " << params.position().transpose());
        ACTS_VERBOSE("  momentum: " << params.momentum().transpose());
        // Push the fitted parameters to the container
        indexedParams.emplace(fitOutput.trackTip, std::move(params));
      } else {
        ACTS_DEBUG("No fitted paramemeters for track " << itrack);
      }
      // Create a PixelMultiTrajectory
      trajectories.emplace_back(std::move(fitOutput.fittedStates),
                                std::move(trackTips), std::move(indexedParams));
    } else {
      ACTS_WARNING("Fit failed for track " << itrack << " with error"
                                           << result.error());
      // Fit failed, but still create an empty PixelMultiTrajectory
      trajectories.push_back(PixelMultiTrajectory());
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return FW::ProcessCode::SUCCESS;
}
