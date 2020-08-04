// This file is part of the Acts project.
//
// Copyright (C) 2020 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cassert>
#include <stdexcept>

#include "ACTFW/Validation/MLTrackClassifier.hpp"
#include "Acts/EventData/MultiTrajectoryHelpers.hpp"

// forward declarations
Acts::ActsVectorXd weightedInput(
    const Acts::ActsMatrixXd& weights, const Acts::ActsVectorXd& input);
Acts::ActsVectorXd reluActivation(const Acts::ActsVectorXd& input);
Acts::ActsVectorXd sigmoidActivation(const Acts::ActsVectorXd& input);

// member functions

FW::MLTrackClassifier::MLTrackClassifier(
    const std::vector<Acts::ActsMatrixXd>& weightsPerLayer)
    : m_weightsPerLayer{weightsPerLayer} {}

FW::MLTrackClassifier::TrackLabels FW::MLTrackClassifier::predictTrackLabel(
    const Acts::MultiTrajectory<SimSourceLink>& multiTraj,
    const size_t& entryIndex, const double& decisionThreshProb) const {
  // check the decision threshold is a probability
  if ((decisionThreshProb < 0.) || (decisionThreshProb > 1.)) {
    throw std::invalid_argument("predictTrackLabel: Decision threshold "
                                "probability is not in [0, 1].");
  }
  // get the trajectory summary info
  auto trajState =
          Acts::MultiTrajectoryHelpers::trajectoryState(multiTraj, entryIndex);
  
  // the vector of input features
  // current neural network was trained on 3 features: hits, outliers, chi2/dof
  Acts::ActsVectorXd inputFeatures(3);
  inputFeatures[0] = trajState.nMeasurements;
  inputFeatures[1] = trajState.nOutliers;
  inputFeatures[2] = trajState.chi2Sum * 1.0 / trajState.NDF;
  
  // linear algebra for layer 1 (hidden layer)

  // weighted input for layer 1
  Acts::ActsVectorXd wInputLayer1 = 
          weightedInput(m_weightsPerLayer[0], inputFeatures);
  // output of layer 1
  Acts::ActsVectorXd outputLayer1 = reluActivation(wInputLayer1);

  // linear algebra for layer 2 (output layer)

  // weighted input for layer 2
  Acts::ActsVectorXd wInputLayer2 = 
          weightedInput(m_weightsPerLayer[1], outputLayer1);
  // output of layer 2
  Acts::ActsVectorXd outputLayer2 = sigmoidActivation(wInputLayer2);

  // the output layer computes how confident the network is that the track is a
  // duplicate; the decision threshold is taken to be the default 50%
  if (outputLayer2[0] > decisionThreshProb) {
      return TrackLabels::duplicate;
  } else {
      return TrackLabels::good;
  }
}

size_t FW::MLTrackClassifier::getNumberofLayers() const {
  return m_weightsPerLayer.size();
}

Acts::ActsMatrixXd FW::MLTrackClassifier::getWeightsAtLayer(
    const size_t& layerIndex) const {
  if (layerIndex >= m_weightsPerLayer.size()) {
    throw std::invalid_argument("getWeightsAtLayer: Layer index invalid.");
  }
  return m_weightsPerLayer[layerIndex];
}

// non-member functions
// these are helper functions internal to predictTrackLabel

// function to compute the weighted input for a layer
Acts::ActsVectorXd weightedInput(
    const Acts::ActsMatrixXd& weights, const Acts::ActsVectorXd& input) {
  // add bias term to input vector
  Acts::ActsVectorXd inputWithBias(input.size() + 1);
  inputWithBias << input, 1.;
  // check that dimensions are compatible for matrix - vector multiplication
  if (weights.cols() != inputWithBias.rows()) {
    throw std::runtime_error("weightedInput: Matrix dimensions incompatible "
                             "with vector dimensions.");
  }
  // weighted_input = weights_matrix * input
  return weights * inputWithBias;
}

// function that applies ReLU activation to the layer input
Acts::ActsVectorXd reluActivation(const Acts::ActsVectorXd& input) {
  // relu(z) = max(0, z)
  Acts::ActsVectorXd reluInput = input.cwiseMax(0);
  return reluInput;
}

// function that applies sigmoid activation to the layer input
Acts::ActsVectorXd sigmoidActivation(const Acts::ActsVectorXd& input) {
  // sigmoid(z) = exp(z) / (1 + exp(z))
  Eigen::ArrayXd expInput = input.array().exp();
  Acts::ActsVectorXd sigmoidInput;
  sigmoidInput = expInput / (expInput + Eigen::ArrayXd::Ones(expInput.size()));
  return sigmoidInput;
}