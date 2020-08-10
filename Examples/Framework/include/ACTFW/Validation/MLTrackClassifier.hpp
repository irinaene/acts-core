// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>

#include "ACTFW/EventData/SimSourceLink.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace FW {

// Classify tracks as good/duplicate/fake using using a deep neural network.
class MLTrackClassifier {
 public:
  /// @brief The labels for track quality
  enum TrackLabels { good, duplicate, fake };

  /// @brief Default constructor
  MLTrackClassifier() = default;

  /// @brief Parametrized constructor
  ///
  /// @param vector that holds the weights matrix for each layer
  MLTrackClassifier(const std::vector<Acts::ActsMatrixXd>& weightsPerLayer);

  /// @brief Default destructor
  ~MLTrackClassifier() = default;

  /// @brief predict the track label
  ///
  /// @param multiTraj The MultiTrajectory object
  /// @param entryIndex The entry index of the trajectory to be classified
  ///
  /// @return the predicted track label of the trajectory
  TrackLabels predictTrackLabel(
      const Acts::MultiTrajectory<SimSourceLink>& multiTraj,
      const size_t& entryIndex, const double& decisionThreshProb) const;

  /// @brief getter for the number of layers in the neural network
  ///
  /// @return the total number of layers in the network
  size_t getNumberofLayers() const;

  /// @brief getter for the weights matrix of a specific layer of the network
  ///
  /// @param layerIndex identifies the location of the layer in the network
  ///
  /// @return matrix of weights specific to the layer indexed by layerIndex
  Acts::ActsMatrixXd getWeightsAtLayer(const size_t& layerIndex) const;

 private:
  /// Vector that holds the weights matrix for each layer, in the order the
  /// layers appear in the network
  std::vector<Acts::ActsMatrixXd> m_weightsPerLayer;
};

}  // namespace FW