// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <limits>
#include <map>
#include <vector>

#include "Acts/Alignment/AlignmentError.hpp"
#include "Acts/Alignment/detail/AlignmentEngine.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/Fitter/detail/KalmanGlobalCovariance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AlignmentDefinitions.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
///
/// @brief Options for align() call
///
/// @tparam fit_options_t The fit options type
///
template <typename fit_options_t>
struct AlignmentOptions {
  /// Deleted default constructor
  AlignmentOptions() = delete;

  /// AlignmentOptions
  ///
  /// @param fOptions The fit options
  /// @param aSurfaces The alignable surfaces
  /// @param chi2CufOff The alignment chi2 tolerance
  /// @param maxIters The alignment maximum iterations
  AlignmentOptions(const fit_options_t& fOptions,
                   const std::vector<const Surface*>& aSurfaces = {},
                   double chi2CutOff = 1, size_t maxIters = 1)
      : fitOptions(fOptions),
        alignableSurfaces(aSurfaces),
        deltaChi2CutOff(chi2CutOff),
        maxIterations(maxIters) {}

  // The fit options
  fit_options_t fitOptions;

  // The surfaces (or detector elements?) to be aligned
  std::vector<const Surface*> alignableSurfaces;

  // The alignment tolerance
  double deltaChi2CutOff = 1;

  // The maximum number of iterations to run alignment
  size_t maxIterations = 5;
};

/// @brief Alignment result struct
///
struct AlignmentResult {
  // The change of alignment parameters
  ActsVectorX<BoundParametersScalar> deltaAlignmentParameters;

  // The change of chi2
  double deltaChi2 = 0;

  // The covariance of alignment parameters
  ActsMatrixX<BoundParametersScalar> alignmentCovariance;

  // The minimized average chi2 per track
  // double averagedChi2 = std::numeric_limits<double>::max();

  // The number of alignment dof
  size_t dof = 0;

  Result<void> result{Result<void>::success()};
};

/// @brief KalmanFitter-based alignment implementation
///
/// @tparam fitter_t Type of the fitter class
template <typename fitter_t>
struct Alignment {
  /// Default constructor is deleted
  Alignment() = delete;

  /// Constructor from arguments
  Alignment(fitter_t fitter, std::unique_ptr<const Logger> logger =
                                 getDefaultLogger("Alignment", Logging::INFO))
      : m_fitter(std::move(fitter)), m_logger(logger.release()) {}

  /// @brief evaluate alignment state for a single track
  ///
  /// @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fit_options_t The fit options type
  ///
  /// @param gctx The current geometry context object
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param fitOptions The fit Options steering the fit
  /// @param idxedAlignSurfaces The idxed surfaces to be aligned
  ///
  /// @result The alignment state for a single track
  template <typename source_link_t, typename start_parameters_t,
            typename fit_options_t>
  Result<detail::TrackAlignmentState> evaluateTrackAlignmentState(
      const GeometryContext& gctx,
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters, const fit_options_t& fitOptions,
      const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces)
      const {
    // Perform the fit
    auto fitRes = m_fitter.fit(sourcelinks, sParameters, fitOptions);
    if (not fitRes.ok()) {
      return fitRes.error();
    }
    // The fit results
    const auto& fitOutput = fitRes.value();
    // Calculate the global track parameters covariance with the fitted track
    const auto& globalTrackParamsCov = detail::globalTrackParametersCovariance(
        fitOutput.fittedStates, fitOutput.trackTip);
    // Calculate the alignment state
    const auto alignState = detail::trackAlignmentState(
        gctx, fitOutput.fittedStates, fitOutput.trackTip, globalTrackParamsCov,
        idxedAlignSurfaces);
    if (alignState.alignmentDof == 0) {
      return AlignmentError::NoAlignmentDofOnTrack;
    }
    return alignState;
  }

  /// @brief update the alignment parameters
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param fitOptions The fit Options steering the fit
  /// @param idxedAlignSurfaces The indexed surfaces to be aligned
  /// @param alignResult [in, out] The aligned result
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  Result<void> updateAlignmentParameters(
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const fit_options_t& fitOptions,
      const std::unordered_map<const Surface*, size_t>& idxedAlignSurfaces,
      AlignmentResult& alignResult) const {
    // The number of trajectories must be eual to the number of starting
    // parameters
    assert(trajectoryCollection.size() == startParametersCollection.size());

    // The total alignment degree of freedom
    alignResult.dof = idxedAlignSurfaces.size() * eAlignmentParametersSize;
    // Initialize derivative of chi2 w.r.t. aligment parameters for all tracks
    ActsVectorX<BoundParametersScalar> sumChi2Derivative =
        ActsVectorX<BoundParametersScalar>::Zero(alignResult.dof);
    ActsMatrixX<BoundParametersScalar> sumChi2SecondDerivative =
        ActsMatrixX<BoundParametersScalar>::Zero(alignResult.dof,
                                                 alignResult.dof);
    // Copy the fit options
    fit_options_t fitOptionsWithRefSurface = fitOptions;
    // Calculate contribution to chi2 derivatives from all input trajectories
    for (unsigned int iTraj = 0; iTraj < trajectoryCollection.size(); iTraj++) {
      const auto& sourcelinks = trajectoryCollection.at(iTraj);
      const auto& sParameters = startParametersCollection.at(iTraj);
      // Set the target surface
      fitOptionsWithRefSurface.referenceSurface =
          &sParameters.referenceSurface();
      // The result for one single track
      auto evaluateRes = evaluateTrackAlignmentState(
          fitOptions.geoContext, sourcelinks, sParameters,
          fitOptionsWithRefSurface, idxedAlignSurfaces);
      if (evaluateRes.ok()) {
        const auto& alignState = evaluateRes.value();
        for (const auto& [rowSurface, rows] : alignState.alignedSurfaces) {
          const auto& [dstRow, srcRow] = rows;
          // Fill the results into full chi2 derivative matrixs
          sumChi2Derivative.segment<eAlignmentParametersSize>(
              dstRow * eAlignmentParametersSize) +=
              alignState.alignmentToChi2Derivative.segment(
                  srcRow * eAlignmentParametersSize, eAlignmentParametersSize);

          for (const auto& [colSurface, cols] : alignState.alignedSurfaces) {
            const auto& [dstCol, srcCol] = cols;
            sumChi2SecondDerivative
                .block<eAlignmentParametersSize, eAlignmentParametersSize>(
                    dstRow * eAlignmentParametersSize,
                    dstCol * eAlignmentParametersSize) +=
                alignState.alignmentToChi2SecondDerivative.block(
                    srcRow * eAlignmentParametersSize,
                    srcCol * eAlignmentParametersSize, eAlignmentParametersSize,
                    eAlignmentParametersSize);
          }
        }
      }
    }
    // Check if the chi2 second derivative matrix inverse is valid
    if (sumChi2SecondDerivative.inverse().hasNaN()) {
      ACTS_ERROR("Chi2 second derivative has NaN");
      return AlignmentError::AlignmentParametersUpdateFailure;
    }

    // Initialize the alignment results
    alignResult.deltaAlignmentParameters =
        ActsVectorX<BoundParametersScalar>::Zero(alignResult.dof);
    alignResult.alignmentCovariance = ActsMatrixX<BoundParametersScalar>::Zero(
        alignResult.dof, alignResult.dof);
    // Alignment parameters change
    alignResult.deltaAlignmentParameters =
        sumChi2SecondDerivative.inverse() * sumChi2Derivative;
    // Alignment parameters covariance
    alignResult.alignmentCovariance = 2 * sumChi2SecondDerivative.inverse();
    // chi2 change
    alignResult.deltaChi2 = 0.5 * sumChi2Derivative.transpose() *
                            alignResult.deltaAlignmentParameters;

    // @Todo: the change of alignment parameters should be passed back to the
    // aligned surfaces

    return Result<void>::success();
  }

  /// @brief Alignment implementation
  ///
  /// @tparam trajectory_container_t The trajectories container type
  /// @tparam start_parameters_t The initial parameters container type
  /// @tparam fit_options_t The fit options type
  ///
  /// @param trajectoryCollection The collection of trajectories as input of
  /// fitting
  /// @param startParametersCollection The collection of starting parameters as
  /// input of fitting
  /// @param alignOptions The alignment options
  ///
  /// @result The alignment result
  template <typename trajectory_container_t,
            typename start_parameters_container_t, typename fit_options_t>
  Result<AlignmentResult> align(
      const trajectory_container_t& trajectoryCollection,
      const start_parameters_container_t& startParametersCollection,
      const AlignmentOptions<fit_options_t>& alignOptions) const {
    // Construct an AlignmentResult object
    AlignmentResult alignRes;

    // Assign index to the alignable surface
    std::unordered_map<const Surface*, size_t> idxedAlignSurfaces;
    for (unsigned int iSurface = 0;
         iSurface < alignOptions.alignableSurfaces.size(); iSurface++) {
      idxedAlignSurfaces.emplace(alignOptions.alignableSurfaces.at(iSurface),
                                 iSurface);
    }

    // Start the iteration to minimize the chi2
    bool converged = false;
    for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
      // Perform the fit to the trajectories and update alignment parameters
      auto updateRes = updateAlignmentParameters(
          trajectoryCollection, startParametersCollection,
          alignOptions.fitOptions, idxedAlignSurfaces, alignRes);
      if (not updateRes.ok()) {
        ACTS_ERROR("Update alignment parameters failed: " << updateRes.error());
        return updateRes.error();
      }
      ACTS_VERBOSE("iIter = " << iIter
                              << ", deltaChi2 = " << alignRes.deltaChi2);
      // Check if it has converged against the provided precision
      if (alignRes.deltaChi2 <= alignOptions.deltaChi2CutOff) {
        converged = true;
        break;
      }
    }
    // Alignment failure if not converged
    if (not converged) {
      ACTS_ERROR("Alignment is not converged.");
      alignRes.result = AlignmentError::ConvergeFailure;
    }

    return alignRes;
  }

 private:
  // The fitter
  fitter_t m_fitter;

  /// Logger getter to support macros
  const Logger& logger() const { return *m_logger; }

  /// Owned logging instance
  std::shared_ptr<const Logger> m_logger;
};
}  // namespace Acts
