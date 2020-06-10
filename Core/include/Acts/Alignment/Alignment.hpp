// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <limits>
#include <map>
#include <vector>

#include "Acts/Alignment/AlignmentError.hpp"
#include "Acts/Alignment/detail/AlignmentEngine.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

namespace Acts {
///
/// @brief Options for alignment() call
///
struct AlignmentOptions {
  // The surfaces (or detector elements?) to be aligned
  std::vector<const Surface*> alignableSurfaces;

  // The alignment tolerance
  double deltaChi2Cutoff = 1;

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
  /// @tparam fitter_options_t Type of the kalman fitter options
  ///
  /// @param sourcelinks The fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param fitOptions KalmanOptions steering the fit
  /// @param aSurfaces The indexed surfaces to be aligned
  ///
  /// @param result The alignment state for a single track
  template <typename source_link_t, typename start_parameters_t,
            typename fitter_options_t>
  Result<TrackAlignmentState> evaluateTrackAlignmentState(
      const std::vector<source_link_t>& sourcelinks,
      const start_parameters_t& sParameters, const fitter_options_t& fitOptions,
      const std::unordered_map<const Surface*, size_t>& aSurfaces) const {
    // Perform the fit
    auto fitRes = m_fitter.fit(sourcelinks, sParameters, fitOptions);
    if (not fitRes.ok()) {
      return fitRes.error();
    }
    // The fit results
    const auto& fitOutput = fitRes.value();
    // Calculate the global track parameters covariance with the fitted track
    const auto& globalTrackParamsCov = detail::globalTrackParametersCov(
        fitOutput.fittedStates, fitOutput.trackTip);
    // Calculate the alignment state
    const auto alignState =
        detail::trackAlignmentState(fitOutput.fittedStates, fitOutput.trackTip,
                                    globalTrackParamsCov, aSurfaces);
    if (alignState.alignmentDof == 0) {
      return AlignmentError : NoAlignmentDofOnTrack;
    }
    return alignState;
  }

  /// @brief update the alignment parameters
  ///
  // @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fitter_options_t Type of the kalman fitter options
  ///
  /// @param inputs The pair of input source links and initial track parameters
  /// used to run fitting for one trajectory
  /// @param fitOptions KalmanOptions steering the fit
  /// @param aSurfaces The indexed surfaces to be aligned
  /// @param alignResult [in, out] The aligned result
  template <typename source_link_t, typename start_parameters_t,
            typename fitter_options_t>
  Result<void> updateAlignmentParameters(
      const std::vector<std::pair<const std::vector<source_link_t>,
                                  const start_parameters_t>>& inputs,
      const fitter_options_t& fitOptions,
      const std::unordered_map<const Surface*, size_t>& aSurfaces,
      AlignmentResult& alignResult) const {
    // The total alignment degree of freedom
    alignResult.dof = aSurfaces.size() * eAlignmentParametersSize;
    // Initialize derivative of chi2 w.r.t. aligment parameters for all tracks
    ActsVectorX<BoundParametersScalar> sumChi2Derivative =
        ActsVectorX<BoundParametersScalar>::Zero(alignResult.dof);
    ActsMatrixX<BoundParametersScalar> sumChi2SecondDerivative =
        ActsMatrixX<BoundParametersScalar>::Zero(alignResult.dof,
                                                 alignResult.dof);
    // Calculate contribution to chi2 derivatives from all input trajectories
    for (const auto& [sourcelinks, sParameters] : inputs) {
      // The result for one single track
      auto eRes = evaluateTrackAlignmentState(sourcelinks, sParameters,
                                              fitOptions, aSurfaces);
      if (eRes.ok()) {
        const auto& alignState = eRes.value();
        for (const auto& [rowSurface, [ dstRow, srcRow ]] :
             alignState.alignedSurfaces) {
          // Fill the results into full chi2 derivative matrixs
          sumChi2Derivative.segment<eAlignmentParameters>(
              dstRow * eAlignmentParameters) +=
              alignState.alignmentToChi2Derivative
                  .segment<eAlignmentParameters>(srcRow * eAlignmentParameters);

          for (const auto& [colSurface, [ dstCol, srcCol ]] :
               alignState.alignedSurfaces) {
            sumChi2SecondDerivative
                .block<eAlignmentParameters, eAlignmentParameters>(
                    dstRow * eAlignmentParameters,
                    dstCol * eAlignmentParameters) +=
                alignState.alignmentToChi2SecondDerivative
                    .block<eAlignmentParameters, eAlignmentParameters>(
                        srcRow * eAlignmentParameters,
                        srcCol * eAlignmentParameters);
          }
        }
      }
    }
    // Check if the chi2 second derivative matrix inverse is valid
    if (sumChi2SecondDerivative.inverse().hasNaN()) {
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
    alignResult.deltaChi2 =
        0.5 * sumChi2Derivative.transpose() * deltaAlignmentParameters;

    // @Todo: the change of alignment parameters should be passed back to the
    // aligned surfaces

    return Result<void>::success();
  }

  /// @brief Alignment implementation
  ///
  // @tparam source_link_t Source link type identifying uncalibrated input
  /// measurements.
  /// @tparam start_parameters_t Type of the initial parameters
  /// @tparam fitter_options_t Type of the kalman fitter options
  ///
  /// @param inputs The pair of input source links and initial track parameters
  /// used to run fitting for one trajectory
  /// @param fitOptions KalmanOptions steering the fit
  /// @param alignOptions AlignmentOptions steering the alignment
  ///
  /// @param result The alignment result
  template <typename source_link_t, typename start_parameters_t,
            typename fitter_options_t>
  Result<AlignmentResult> alignment(
      const std::vector<std::pair<const std::vector<source_link_t>,
                                  const start_parameters_t>>& inputs,
      const fitter_options_t& fitOptions,
      const AlignmentOptions& alignOptions) const {
    // Construct an AlignmentResult object
    AlignmentResult alignRes;

    // Assign index to the alignable surface
    std::unordered_map<const Surface*, size_t> aSurfaces;
    for (unsigned int iSurface = 0;
         iSurface < alignOptions.alignableSurfaces.size(); iSurface++) {
      aSurfaces.emplace(alignOptions.alignableSurfaces.at(iSurface), iSurface);
    }

    // Start the iteration to minimize the chi2
    bool alignSucceed = false;
    for (unsigned int iIter = 0; iIter < alignOptions.maxIterations; iIter++) {
      auto uRes =
          updateAlignmentParameters(inputs, fitOptions, aSurfaces, alignRes);
      if (not uRes.ok()) {
        return uRes.error();
      }
      // Check if it has converged against the provided precision
      if (alignRes.detalChi2 <= alignOptions.deltaChi2CutOff) {
        alignSucceed = true;
        break;
      }
    }
    // Alignment failure if not converged
    if (not alignSucceed) {
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
