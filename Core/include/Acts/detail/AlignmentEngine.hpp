// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
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
/// derivatives matrix
///
/// @tparam source_link_t The source link type of the trajectory
///
/// @param gctx The current geometry context object
/// @param multiTraj The MultiTrajectory containing the trajectory to be
/// investigated
/// @param entryIndex The trajectory entry index
/// @param fullTrackParamsCov The global track parameters covariance for a
/// single track. This contains all track states including those non-measurment
/// states. So re-fill to shrink the matrix could be needed.
/// @param alignSurfaces The indexed surfaces to be aligned
///
/// @return The first and second derivative of chi2 w.r.t. alignment parameters
template <typename source_link_t>
std::tuple<ActsVectorX<BoundParametersScalar>,
           ActsMatrixX<BoundParametersScalar>>
singleTrackAlignmentToChi2Derivatives(
    const GeometryContext& gctx,
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex,
    const ActsMatrixX<BoundParametersScalar>& fullTrackParamsCov,
    const std::unordered_map<const Surface*, size_t>& alignSurfaces) {
  // Remember the indices of to-be-aligned surfaces relevant with this track
  std::vector<size_t> surfaceIndices;
  surfaceIndices.reserve(15);
  // The storage index of the relevant measurement state
  std::vector<size_t> measurementIndices;
  measurementIndices.reserve(15);
  // The ismoothed index of the relevant measurement state (to help retrieve
  // elements stored in global track parameters covariance)
  std::vector<size_t> iSmoothedIndices;
  iSmoothedIndices.reserve(15);
  // The dimension of provided global track parameters covariance should be same
  // as eBoundParametersSize * nSmoothedStates
  size_t nSmoothedStates = 0;
  size_t totalMeasurementDim = 0;
  // Visit the track states on the trajectory
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    if (state.hasSmoothed()) {
      nSmoothedStates++;
    }
    // Only measurement states matter
    if (not state.typeFlags().test(TrackStateFlag::MeasurementFlag)) {
      return true;
    }
    // Check if the reference surface is to be aligned
    const auto surface = &state.referenceSurface();
    auto it = alignSurfaces.find(surface);
    if (it == alignSurfaces.end()) {
      return true;
    }
    // @Todo:add struct for those info
    surfaceIndices.push_back(it->second);
    measurementIndices.push_back(state.index());
    iSmoothedIndices.push_back(nSmoothedStates);
    totalMeasurementDim += state.calibratedSize();
    return true;
  });

  // @Todo: put the following into separate functions
  // Check if the size of provided global track parameters covariance is correct
  assert(trackParamsCov.rows() == trackParamsCov.cols() and
         trackParamsCov.rows() == eBoundParametersSize * nSmoothedStates);

  // Prepare the alignment matrixs which are composed by components from the
  // selected track states
  size_t alignDof = eAlignmentParametersSize * surfaceIndices.size();
  size_t totalParametersDim = eBoundParametersSize * measurementIndices.size();
  // The measurement covariance
  ActsMatrixX<BoundParametersScalar> measurementCovariance;
  measurementCovariance.resize(totalMeasurementDim, totalMeasurementDim);
  // The bound parameters to measurement projection matrix
  ActsMatrixX<BoundParametersScalar> projectionMatrix;
  projectionMatrix.resize(totalMeasurementDim, totalParametersDim);
  // The derivative of residual w.r.t. alignment parameters
  ActsMatrixX<BoundParametersScalar> alignmentToResidualDerivative;
  alignmentToResidualDerivative.resize(totalMeasurementDim, alignDof);
  // The track parameters covariance
  ActsMatrixX<BoundParametersScalar> trackParametersCovariance;
  trackParametersCovariance.resize(totalParametersDim, totalParametersDim);
  // The residual
  ActsVectorX<ParValue_t> residual;
  residual.resize(totalMeasurementDim);

  size_t iMeasurement = totalMeasurementDim;
  size_t iParams = totalParametersDim;
  for (size_t iState = 0; iState < measurementIndices.size(); iState++) {
    const auto& state = multiTraj.getTrackState(measurementIndices.at(iState));
    size_t measdim = state.calibratedSize();
    iMeasurement -= measdim;
    iParams -= eBoundParametersSize;
    // (a) Get and fill the measurement covariance matrix
    ActsSymMatrixD<measdim> thisMeasCovariance =
        state.calibratedCovariance().template topLeftCorner<measdim, measdim>();
    measurementCovariance.block<measdim, measdim>(iMeasurement, iMeasurement) =
        thisMeasCovariance;

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

    // (e) Shrink the provided full track parameters covariance matrix
    size_t iSmoothed = iSmoothedIndices.at(iState);
    for (size_t iPrevState = 0; iPrevState < iSmoothedIndices.size();
         iPrevState++) {
      size_t iPrev = iSmoothedIndices.at(iPrevState);
      size_t iRow = (nSmoothedStates - iSmoothed) * eBoundParametersSize;
      size_t iCol = (nSmoothedStates - iPrev) * eBoundParametersSize;
      BoundSymMatrix thisParamsCovariance =
          fullTrackParametersCov
              .block<eBoundParametersSize, eBoundParametersSize>(iRow, iCol);
      trackParametersCovariance
          .block<eBoundParametersSize, eBoundParametersSize>(
              iParams, iPrevState * eBoundParametersSize) =
          thisParamsCovariance;
    }
  }

  // Calculate the chi2 derivatives
  ActsVectorX<BoundParametersScalar> chi2FirstDerivative;
  ActsMatrixX<BoundParametersScalar> chi2SecondDerivative;
  chi2FirstDerivative.resize(alignDof);
  chi2SecondDerivative.resize(alignDof, alignDof);
  chi2firstDerivative =
      2 * alignmentToResidualDerivative.transpose() *
      measurementCovariance.inverse() *
      (measurementCovariance - projectionMatrix * trackParametersCovariance *
                                   projectionMatrix.transpose()) *
      measurementCovariance.inverse() * residual;
  chi2SecondDerivative =
      2 * alignmentToResidualDerivative.transpose() *
      measurementCovariance.inverse() *
      (measurementCovariance - projectionMatrix * trackParametersCovariance *
                                   projectionMatrix.transpose()) *
      measurementCovariance.inverse() * alignmentToResidualDerivative;

  return std::make_pair(chi2FirstDerivative, chi2SecondDerivative);
}

}  // namespace detail
}  // namespace Acts
