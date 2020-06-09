// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace Acts {
namespace detail {
/// Calculate the first and second derivative of chi2 w.r.t. alignment
/// parameters for a single track
///
/// Suppose there are n measurements on the track, and m (m<=n) of them are on
/// to-be-aligned surface, then eAlignmentParametersSize*m alignment parameters
/// will be involved for this particular track, i.e. this track will contribute
/// to at most eAlignmentParametersSize*m*2 elements of the full chi2
/// second derivative matrix
///
/// @tparam source_link_t The source link type of the trajectory
/// @tparam parameters_t The track parameters type
///
/// @param multiTraj The MultiTrajectory containing the trajectory to be
/// investigated
/// @param entryIndex The trajectory entry index
/// @param fullTrackParamsCov The global track parameters covariance for a
/// single track. This contains all smoothed track states including those
/// non-measurement states. Selection of certain rows/columns for
/// alignment-relevant states is needed.
/// @param alignSurfaces The indexed surfaces to be aligned
///
/// @return The first and second derivative of chi2 w.r.t. alignment parameters
/// and indices of surfaces for those alignment parameters
template <typename source_link_t, typename parameters_t = BoundParameters>
std::tuple<ActsVectorX<BoundParametersScalar>,
           ActsMatrixX<BoundParametersScalar>, std::vector<size_t>>
singleTrackAlignmentToChi2Derivatives(
    const MultiTrajectory<source_link_t>& multiTraj, const size_t& entryIndex,
    const ActsMatrixX<BoundParametersScalar>& fullTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& alignSurfaces) {
  using CovMatrix_t = typename parameters_t::CovMatrix_t;

  // Remember the indices of to-be-aligned surfaces relevant with this track
  std::vector<size_t> surfaceIndices;
  surfaceIndices.reserve(15);
  // Remember the index within the trajectory and smoothing order of the
  // alignment-relevant measurement states (The smoothing order is needed to
  // retrieve relevant elements stored in the full global track parameters
  // covariance)
  std::vector<std::pair<size_t, size_t>> alignStates;
  alignStates.reserve(15);
  // Number of smoothed states on the trajectory
  size_t nSmoothedStates = 0;
  // Dimension of alignment-relevant measurements on the trajectory
  size_t alignMeasurementsDim = 0;

  // Visit the track states on the trajectory
  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    // Remember the number of smoothed states
    if (ts.hasSmoothed()) {
      nSmoothedStates++;
    }
    // Only measurement states matter
    if (not ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
      return true;
    }
    // Check if the reference surface is to be aligned
    const auto surface = &ts.referenceSurface();
    auto it = alignSurfaces.find(surface);
    if (it == alignSurfaces.end()) {
      return true;
    }
    // Remember the index of this surface within the to-be-aligned surfaces pool
    surfaceIndices.push_back(it->second);
    // Rember the index of the state within the trajectory and its smoothing
    // order
    alignStates.push_back({ts.index(), nSmoothedStates});
    // Add up measurement dimension
    alignMeasurementsDim += ts.calibratedSize();
    return true;
  });

  // The dimension of provided global track parameters covariance should be same
  // as eBoundParametersSize * nSmoothedStates
  assert(trackParamsCov.rows() == trackParamsCov.cols() and
         trackParamsCov.rows() == eBoundParametersSize * nSmoothedStates);

  // @Todo: put the following into separate functions
  // Prepare the alignment matrixs composed by components from the
  // alignment-relevant track states
  // The alignment degree of freedom
  size_t alignDof = eAlignmentParametersSize * surfaceIndices.size();
  // Dimension of alignment-relevant parameters on the trajectory
  size_t alignParametersDim = eBoundParametersSize * alignStates.size();

  // The measurement covariance
  ActsMatrixX<BoundParametersScalar> measurementCovariance;
  measurementCovariance.resize(alignMeasurementsDim, alignMeasurementsDim);
  // The bound parameters to measurement projection matrix
  ActsMatrixX<BoundParametersScalar> projectionMatrix;
  projectionMatrix.resize(alignMeasurementsDim, alignParametersDim);
  // The derivative of residual w.r.t. alignment parameters
  ActsMatrixX<BoundParametersScalar> alignmentToResidualDerivative;
  alignmentToResidualDerivative.resize(alignMeasurementsDim, alignDof);
  // The track parameters covariance
  ActsMatrixX<BoundParametersScalar> trackParametersCovariance;
  trackParametersCovariance.resize(alignParametersDim, alignParametersDim);
  // The residual
  ActsVectorX<ParValue_t> residual;
  residual.resize(alignMeasurementsDim);

  size_t iMeasurement = alignMeasurementsDim;
  size_t iParams = alignParametersDim;
  // Loop over the alignment-relevant states to fill those alignment matrixs
  for (size_t iRowState = 0; iRowState < alignStates.size(); iRowState++) {
    const auto& [index, iRowSmoothed] = alignStates.at(iRowState);
    const auto& state = multiTraj.getTrackState(index);
    size_t measdim = state.calibratedSize();
    // The index of current measurement and parameter
    iMeasurement -= measdim;
    iParams -= eBoundParametersSize;
    // (a) Get and fill the measurement covariance matrix
    ActsSymMatrixD<measdim> measCovariance =
        state.calibratedCovariance().template topLeftCorner<measdim, measdim>();
    measurementCovariance.block<measdim, measdim>(iMeasurement, iMeasurement) =
        measCovariance;

    // (b) Get and fill the bound parameters to measurement projection matrix
    const ActsMatrixD<measdim, eBoundParametersSize> H =
        state.projector()
            .template topLeftCorner<measdim, eBoundParametersSize>();
    projectionMatrix.block<measdim, eBoundParametersSize>(iMeasurement,
                                                          iParams) = H;

    // (c) Get and fill the residual
    residual.segment<measdim>(iMeasurement) =
        state.calibrated().template head<meas>() - H * state.filtered();

    // (d) @Todo: Get the derivative of alignment parameters w.r.t. measurement
    // or residual

    // (e) Get and fill the track parameters covariance matrix (for only
    // alignment-relevant states)
    // @Todo: add helper function to select rows/columns of a matrix
    for (size_t iColState = 0; iColState < alignStates.size(); iColState++) {
      size_t iColSmoothed = alignStates.at(iColState).second;
      size_t iRowFull = (nSmoothedStates - iRowSmoothed) * eBoundParametersSize;
      size_t iColFull = (nSmoothedStates - iColSmoothed) * eBoundParametersSize;
      // Retrieve the block from the provided full covariance matrix
      CovMatrix_t correlation =
          fullTrackParamsCov.block<eBoundParametersSize, eBoundParametersSize>(
              iRowFull, iColFull);
      // Fill the block of the target covariance matrix
      trackParametersCovariance
          .block<eBoundParametersSize, eBoundParametersSize>(
              iParams,
              alignParametersDim - (iColState + 1) * eBoundParametersSize) =
          correlation;
    }
  }

  // Calculate the chi2 derivatives based on the alignment matrixs
  ActsVectorX<BoundParametersScalar> chi2FirstDerivative;
  ActsMatrixX<BoundParametersScalar> chi2SecondDerivative;
  chi2FirstDerivative.resize(alignDof);
  chi2SecondDerivative.resize(alignDof, alignDof);
  // The covariance of residual
  ActsMatrixX<BoundParametersScalar> residualCovariance;
  residualCovariance.resize(alignMeasurementsDim, alignMeasurementsDim);
  residualCovariance = measurementCovariance - projectionMatrix *
                                                   trackParametersCovariance *
                                                   projectionMatrix.transpose();

  chi2firstDerivative = 2 * alignmentToResidualDerivative.transpose() *
                        measurementCovariance.inverse() * residualCovariance *
                        measurementCovariance.inverse() * residual;
  chi2SecondDerivative = 2 * alignmentToResidualDerivative.transpose() *
                         measurementCovariance.inverse() * residualCovariance *
                         measurementCovariance.inverse() *
                         alignmentToResidualDerivative;

  return std::make_tuple(chi2FirstDerivative, chi2SecondDerivative,
                         surfaceIndices);
}

}  // namespace detail
}  // namespace Acts
