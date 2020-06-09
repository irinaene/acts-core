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

#include <unordered_map>

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
/// @param globalTrackParamsCov The global track parameters covariance for a
/// single track and the starting row/column for smoothed states. This contains
/// all smoothed track states including those non-measurement states. Selection
/// of certain rows/columns for measurement states is needed.
/// @param alignableSurfaces The indexed surfaces to be aligned
///
/// @return The first and second derivative of chi2 w.r.t. alignment parameters
/// and indices of surfaces for those alignment parameters
template <typename source_link_t, typename parameters_t = BoundParameters>
std::tuple<ActsVectorX<BoundParametersScalar>,
           ActsMatrixX<BoundParametersScalar>, ActsVectorX<size_t>>
singleTrackAlignmentToChi2Derivatives(
    const MultiTrajectory<source_link_t>& multiTraj, const size_t& entryIndex,
    const std::pair<ActsMatrixX<BoundParametersScalar>,
                    std::unordered_map<size_t, size_t>>& globalTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& alignableSurfaces) {
  using CovMatrix_t = typename parameters_t::CovMatrix_t;

  // Remember the index within the trajectory and its reference surface index
  // within the to-be-aligned surfaces pool of the measurement states
  std::vector<std::pair<size_t, size_t>> measurementStates;
  measurementStates.reserve(15);
  // Number of smoothed states on the track
  size_t nSmoothedStates = 0;
  // Dimension of measurements on the track
  size_t trackMeasurementsDim = 0;
  // Number of to-be-aligned surfaces on the track
  size_t nAlignSurfaces = 0;

  // Visit the track states on the track
  multiTraj.visitBackwards(entryIndex, [&](const auto& ts) {
    // Remember the number of smoothed states
    if (ts.hasSmoothed()) {
      nSmoothedStates++;
    }
    // Only measurement states matter (we can't align non-measurement states,
    // no?)
    if (not ts.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
      return true;
    }
    // Check if the reference surface is to be aligned
    // @Todo: consider the case when some of the Dofs are fixed for one surface
    size_t surfaceIndex = SIZE_MAX;
    const auto surface = &ts.referenceSurface();
    auto it = alignableSurfaces.find(surface);
    if (it != alignableSurfaces.end()) {
      surfaceIndex = it->second;
      nAlignSurfaces++;
    }
    // Rember the index of the state within the trajectory and the index of the
    // surface within the to-be-aligned surfaces pool
    measurementStates.push_back({ts.index(), surfaceIndex});
    // Add up measurement dimension
    trackMeasurementsDim += ts.calibratedSize();
    return true;
  });

  // The dimension of provided global track parameters covariance should be same
  // as eBoundParametersSize * nSmoothedStates
  assert(trackParamsCov.rows() == trackParamsCov.cols() and
         trackParamsCov.rows() == eBoundParametersSize * nSmoothedStates);

  // @Todo: put the following into separate functions
  // Prepare the alignment matrixs composed by components from the measurement
  // states The alignment degree of freedom
  size_t alignDof = eAlignmentParametersSize * nAlignSurfaces;
  // Dimension of global track parameters (from only measurement states) on the
  // track
  size_t trackParametersDim = eBoundParametersSize * measurementStates.size();

  // The indices of aligned surfaces on the track
  ActsVectorX<size_t> surfaceIndices =
      ActsVectorX<size_t>::Zero(nAlignSurfaces);
  // The measurement covariance
  ActsMatrixX<BoundParametersScalar> measurementCovariance =
      ActsMatrixX<BoundParametersScalar>::Zero(trackMeasurementsDim,
                                               trackMeasurementsDim);
  // The bound parameters to measurement projection matrix
  ActsMatrixX<BoundParametersScalar> projectionMatrix =
      ActsMatrixX<BoundParametersScalar>::Zero(trackMeasurementsDim,
                                               trackParametersDim);
  // The derivative of residual w.r.t. alignment parameters
  ActsMatrixX<BoundParametersScalar> alignmentToResidualDerivative =
      ActsMatrixX<BoundParametersScalar>::Zero(trackMeasurementsDim, alignDof);
  // The track parameters covariance
  ActsMatrixX<BoundParametersScalar> trackParametersCovariance =
      ActsMatrixX<BoundParametersScalar>::Zero(trackParametersDim,
                                               trackParametersDim);
  // The residual
  ActsVectorX<ParValue_t> residual =
      ActsVectorX<ParValue_t>::Zero(trackMeasurementsDim);

  // Unpack global track parameters covariance and the starting row/column for
  // all smoothed states
  const auto& [sourceTrackParamsCov, stateRowIndices] = globalTrackParamsCov;
  // Loop over the measurement states to fill those alignment matrixs
  // This is done in reverse order
  size_t iMeasurement = trackMeasurementsDim;
  size_t iParams = trackParametersDim;
  size_t iSurface = nAlignSurfaces;
  for (const auto& [rowStateIndex, surfaceIndex] : measurementStates) {
    const auto& state = multiTraj.getTrackState(rowStateIndex);
    size_t measdim = state.calibratedSize();
    // Update index of current measurement and parameter
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
    if (surfaceIndex != SIZE_MAX) {
      iSurface -= 1;
      surfaceIndices.row(iSurface) = surfaceIndex;
    }

    // (e) Extract and fill the track parameters covariance matrix for only
    // measurement states
    // @Todo: add helper function to select rows/columns of a matrix
    for (size_t iColState = 0; iColState < measurementStates.size();
         iColState++) {
      size_t colStateIndex = measurementStates.at(iColState).first;
      // Retrieve the block from the source covariance matrix
      CovMatrix_t correlation =
          sourceTrackParamsCov
              .block<eBoundParametersSize, eBoundParametersSize>(
                  stateRowIndices.at(rowStateIndex),
                  stateRowIndices.at(colStateIndex));
      // Fill the block of the target covariance matrix
      size_t iCol = trackParametersDim - (iColState + 1) * eBoundParametersSize;
      trackParametersCovariance
          .block<eBoundParametersSize, eBoundParametersSize>(iParams, iCol) =
          correlation;
    }
  }

  // Calculate the chi2 derivatives based on the alignment matrixs
  ActsVectorX<BoundParametersScalar> chi2FirstDerivative =
      ActsVectorX<BoundParametersScalar>::Zero(alignDof);
  ActsMatrixX<BoundParametersScalar> chi2SecondDerivative =
      ActsMatrixX<BoundParametersScalar>::Zero(alignDof, alignDof);
  // The covariance of residual
  ActsMatrixX<BoundParametersScalar> residualCovariance =
      ActsMatrixX<BoundParametersScalar>::Zero(trackMeasurementsDim,
                                               trackMeasurementsDim);
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
