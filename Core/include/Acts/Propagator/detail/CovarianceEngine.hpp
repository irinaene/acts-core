// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <functional>
#include <variant>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Propagator/StepperState.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @brief These functions perform the transport of a covariance matrix using
/// given Jacobians. The required data is provided by a @p StepperState object
/// with some additional data. Since this is a purely algebraic problem the
/// calculations are identical for @c StraightLineStepper and @c EigenStepper.
/// As a consequence the methods can be located in a seperate file.
namespace detail {

using Jacobian = std::variant<BoundMatrix, FreeToBoundMatrix, BoundToFreeMatrix, FreeMatrix>;

/// Create and return the bound state at the current position
///
/// @brief It does not check if the transported state is at the surface, this
/// needs to be guaranteed by the propagator
///
/// @param [in] state State that will be presented as @c BoundState
/// @param [in] surface The surface to which we bind the state
///
/// @return A bound state:
///   - the parameters at the surface
///   - the stepwise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
std::tuple<BoundParameters, Jacobian, double> boundState(
    StepperState& state, const Surface& surface);

/// Create and return a curvilinear state at the current position
///
/// @brief This creates a curvilinear state.
///
/// @param [in] state State that will be presented as @c CurvilinearState
///
/// @return A curvilinear state:
///   - the curvilinear parameters at given position
///   - the stepweise jacobian towards it (from last bound)
///   - and the path length (from start - for ordering)
std::tuple<CurvilinearParameters, Jacobian, double> curvilinearState(
    StepperState& state);

  /// Create and return a free state at the current position
  ///
  /// @brief This creates a free state.
  ///
    /// @tparam result_t Defines the return type
  /// @param [in] state State that will be presented as @c FreeState
  /// @param [in] reinitialize Boolean flag whether reinitialization is needed,
  /// i.e. if this is an intermediate state of a larger propagation
  ///
  /// @return A free state:
  ///   - the free parameters at given position
  ///   - the stepweise jacobian towards it (from last location)
  ///   - and the path length (from start - for ordering)
  std::tuple<FreeParameters, Jacobian, double>
freeState(StepperState& state, bool reinitialize);

/// @brief Method for on-demand transport of the covariance to a new frame at current position in parameter space. It treats different scenarios:
/// - from local to local
/// - from local to global
/// - from global to local
/// - from global to global
///
/// @param [in,out] state The stepper state
/// @param [in] toLocal Boolean flag whether the result is in local parameters
/// @param [in] surface is the surface to which the covariance is
///        forwarded to
/// @note No check is done if the position is actually on the surface
///
/// @return Projection jacobian from global to bound parameters
void covarianceTransport(StepperState& state, bool toLocal = true, const Surface* surface = nullptr);

}  // namespace detail
}  // namespace Acts