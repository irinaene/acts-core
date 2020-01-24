// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/detail/CovarianceEngine.hpp"

namespace Acts {
namespace {
/// Some type defs
using Covariance = std::variant<BoundSymMatrix, FreeSymMatrix>;
using BoundState = std::tuple<BoundParameters, detail::Jacobian, double>;
using CurvilinearState =
    std::tuple<CurvilinearParameters, detail::Jacobian, double>;
using FreeState = std::tuple<FreeParameters, detail::Jacobian, double>;

/// @brief Evaluate the projection Jacobian from free to curvilinear parameters
///
/// @param [in] state State that will be projected
///
/// @return Projection Jacobian
FreeToBoundMatrix freeToCurvilinearJacobian(const StepperState& state) {
  // Optimized trigonometry on the propagation direction
  const double x = state.dir(0);  // == cos(phi) * sin(theta)
  const double y = state.dir(1);  // == sin(phi) * sin(theta)
  const double z = state.dir(2);  // == cos(theta)
  // can be turned into cosine/sine
  const double cosTheta = z;
  const double sinTheta = sqrt(x * x + y * y);
  const double invSinTheta = 1. / sinTheta;
  const double cosPhi = x * invSinTheta;
  const double sinPhi = y * invSinTheta;
  // prepare the jacobian to curvilinear
  FreeToBoundMatrix jacToCurv = FreeToBoundMatrix::Zero();
  if (std::abs(cosTheta) < s_curvilinearProjTolerance) {
    // We normally operate in curvilinear coordinates defined as follows
    jacToCurv(0, 0) = -sinPhi;
    jacToCurv(0, 1) = cosPhi;
    jacToCurv(1, 0) = -cosPhi * cosTheta;
    jacToCurv(1, 1) = -sinPhi * cosTheta;
    jacToCurv(1, 2) = sinTheta;
  } else {
    // Under grazing incidence to z, the above coordinate system definition
    // becomes numerically unstable, and we need to switch to another one
    const double c = sqrt(y * y + z * z);
    const double invC = 1. / c;
    jacToCurv(0, 1) = -z * invC;
    jacToCurv(0, 2) = y * invC;
    jacToCurv(1, 0) = c;
    jacToCurv(1, 1) = -x * y * invC;
    jacToCurv(1, 2) = -x * z * invC;
  }
  // Time parameter
  jacToCurv(5, 3) = 1.;
  // Directional and momentum parameters for curvilinear
  jacToCurv(2, 4) = -sinPhi * invSinTheta;
  jacToCurv(2, 5) = cosPhi * invSinTheta;
  jacToCurv(3, 6) = -invSinTheta;
  jacToCurv(4, 7) = 1.;

  return jacToCurv;
}

/// @brief This function treats the modifications of the jacobian related to the
/// projection onto a surface. Since a variation of the start parameters within
/// a given uncertainty would lead to a variation of the end parameters, these
/// need to be propagated onto the target surface. This an approximated approach
/// to treat the (assumed) small change.
///
/// @param [in] state The current state
/// @param [in] surface The surface onto which the projection should be
/// performed
/// @note The parameter @p surface is only required if projected to bound
/// parameters. In the case of curvilinear parameters the geometry and the
/// position is known and the calculation can be simplified
///
/// @return The projection jacobian from global end parameters to its local
/// equivalent
const FreeToBoundMatrix surfaceDerivative(StepperState& state,
                                          const Surface* surface = nullptr) {
  if (state.jacToGlobal.has_value()) {
    // Set the surface projection contributions
    // If no surface is specified it is curvilinear
    if (surface == nullptr) {
      // Transport the covariance
      const ActsRowVectorD<3> normVec(state.dir);
      const BoundRowVector sfactors =
          normVec *
          (*state.jacToGlobal).template topLeftCorner<3, BoundParsDim>();
      *state.jacToGlobal -= state.derivative * sfactors;
      // Since the jacobian to local needs to calculated for the bound
      // parameters here, it is convenient to do the same here
      return freeToCurvilinearJacobian(state);
    }
    // Else it is bound
    else {
      // Initialize the transport final frame jacobian
      FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
      // Initalize the jacobian to local, returns the transposed ref frame
      auto rframeT = surface->initJacobianToLocal(state.geoContext, jacToLocal,
                                                  state.pos, state.dir);
      // Calculate the form factors for the derivatives
      const BoundRowVector sVec =
          surface->derivativeFactors(state.geoContext, state.pos, state.dir,
                                     rframeT, (*state.jacToGlobal));
      *state.jacToGlobal -= state.derivative * sVec;
      // Return the jacobian to local
      return jacToLocal;
    }
  } else {
    // Set the surface projection contributions
    // If no surface is specified it is curvilinear
    if (surface == nullptr) {
      // Transport the covariance
      const ActsRowVectorD<3> normVec(state.dir);
      const FreeRowVector sfactors =
          normVec * state.jacTransport.template topLeftCorner<3, FreeParsDim>();
      // Since the jacobian to local needs to calculated for the bound
      // parameters here, it is convenient to do the same here
      return freeToCurvilinearJacobian(state) *
             (state.jacTransport - state.derivative * sfactors);
    }
    // Else it is bound
    else {
      // Initialize the transport final frame jacobian
      FreeToBoundMatrix jacToLocal = FreeToBoundMatrix::Zero();
      // Initalize the jacobian to local, returns the transposed ref frame
      auto rframeT = surface->initJacobianToLocal(state.geoContext, jacToLocal,
                                                  state.pos, state.dir);
      // Calculate the form factors for the derivatives
      const FreeRowVector sVec = surface->derivativeFactors(
          state.geoContext, state.pos, state.dir, rframeT, state.jacTransport);
      // Return the jacobian to local
      return jacToLocal * (state.jacTransport - state.derivative * sVec);
    }
  }
}

/// @brief This function reinitialises the @p state member @p jacToGlobal.
///
/// @param [in, out] state The state object
/// @param [in] surface The surface the represents the local parametrisation
/// @note The surface is only required for bound parameters since it serves to
/// derive the jacobian from it. In the case of curvilinear parameters this is
/// not needed and can be evaluated without any surface.
void reinitializeJacToGlobal(StepperState& state,
                             const Surface* surface = nullptr) {
  using VectorHelpers::phi;
  using VectorHelpers::theta;

  // Reset the jacobian
  state.jacToGlobal = BoundToFreeMatrix::Zero();

  // Fill the jacobian to global for next transport
  // If treating curvilinear parameters
  if (surface == nullptr) {
    auto& jac = *state.jacToGlobal;
    // TODO: This was calculated before - can it be reused?
    // Optimized trigonometry on the propagation direction
    const double x = state.dir(0);  // == cos(phi) * sin(theta)
    const double y = state.dir(1);  // == sin(phi) * sin(theta)
    const double z = state.dir(2);  // == cos(theta)
    // can be turned into cosine/sine
    const double cosTheta = z;
    const double sinTheta = sqrt(x * x + y * y);
    const double invSinTheta = 1. / sinTheta;
    const double cosPhi = x * invSinTheta;
    const double sinPhi = y * invSinTheta;

    jac(0, eLOC_0) = -sinPhi;
    jac(0, eLOC_1) = -cosPhi * cosTheta;
    jac(1, eLOC_0) = cosPhi;
    jac(1, eLOC_1) = -sinPhi * cosTheta;
    jac(2, eLOC_1) = sinTheta;
    jac(3, eT) = 1;
    jac(4, ePHI) = -sinTheta * sinPhi;
    jac(4, eTHETA) = cosTheta * cosPhi;
    jac(5, ePHI) = sinTheta * cosPhi;
    jac(5, eTHETA) = cosTheta * sinPhi;
    jac(6, eTHETA) = -sinTheta;
    jac(7, eQOP) = 1;
  }
  // If treating bound parameters
  else {
    Vector2D loc{0., 0.};
    surface->globalToLocal(state.geoContext, state.pos, state.dir, loc);
    BoundVector pars;
    pars << loc[eLOC_0], loc[eLOC_1], phi(state.dir), theta(state.dir),
        state.q / state.p, state.t;
    surface->initJacobianToGlobal(state.geoContext, *state.jacToGlobal,
                                  state.pos, state.dir, pars);
  }
}
}  // namespace

namespace detail {

BoundState boundState(StepperState& state, const Surface& surface) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    // Initialize the transport final frame jacobian
    covarianceTransport(state, surface);
    cov = std::get<BoundSymMatrix>(state.cov);
  }
  // Create the bound parameters
  BoundParameters parameters(state.geoContext, cov, state.pos,
                             state.p * state.dir, state.q, state.t,
                             surface.getSharedPtr());
  // Create the bound state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

CurvilinearState curvilinearState(StepperState& state) {
  // Transport the covariance to here
  std::optional<BoundSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, true);
    cov = std::get<BoundSymMatrix>(state.cov);
  }
  // Create the curvilinear parameters
  CurvilinearParameters parameters(cov, state.pos, state.p * state.dir, state.q,
                                   state.t);
  // Create the curvilinear state
  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

FreeState freeState(StepperState& state) {
  // Transport the covariance to here
  std::optional<FreeSymMatrix> cov = std::nullopt;
  if (state.covTransport) {
    covarianceTransport(state, false);
    cov = std::get<FreeSymMatrix>(state.cov);
  }
  // Create the free parameters
  FreeVector pars;
  pars.template head<3>() = state.pos;
  pars(3) = state.t;
  pars.template segment<3>(4) = state.dir;
  pars(7) = (state.q / state.p);
  FreeParameters parameters(std::move(cov), pars);

  return std::make_tuple(std::move(parameters), state.jacobian,
                         state.pathAccumulated);
}

void covarianceTransport(StepperState& state, const Surface& surface) {
  if (state.jacToGlobal.has_value()) {
    state.jacToGlobal = state.jacTransport * (*state.jacToGlobal);
    state.jacobian =
        BoundMatrix(surfaceDerivative(state, &surface) * (*state.jacToGlobal));
    state.cov =
        BoundSymMatrix(std::get<BoundMatrix>(state.jacobian) *
                       std::get<BoundSymMatrix>(state.cov) *
                       std::get<BoundMatrix>(state.jacobian).transpose());
  } else {
    state.jacobian = FreeToBoundMatrix(surfaceDerivative(state, &surface) *
                                       state.jacTransport);
    state.cov =
        BoundSymMatrix(std::get<FreeToBoundMatrix>(state.jacobian) *
                       std::get<FreeSymMatrix>(state.cov) *
                       std::get<FreeToBoundMatrix>(state.jacobian).transpose());
  }

  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
  reinitializeJacToGlobal(state, &surface);
}

void covarianceTransport(StepperState& state, bool toLocal) {
  if (state.jacToGlobal.has_value()) {
    state.jacToGlobal = state.jacTransport * (*state.jacToGlobal);
    if (toLocal) {
      state.jacobian =
          BoundMatrix(surfaceDerivative(state) * (*state.jacToGlobal));
      state.cov =
          BoundSymMatrix(std::get<BoundMatrix>(state.jacobian) *
                         std::get<BoundSymMatrix>(state.cov) *
                         std::get<BoundMatrix>(state.jacobian).transpose());
    } else {
      state.jacobian = BoundToFreeMatrix(*state.jacToGlobal);
      state.cov = FreeSymMatrix(
          std::get<BoundToFreeMatrix>(state.jacobian) *
          std::get<BoundSymMatrix>(state.cov) *
          std::get<BoundToFreeMatrix>(state.jacobian).transpose());
    }
  } else {
    if (toLocal) {
      state.jacobian =
          FreeToBoundMatrix(surfaceDerivative(state) * state.jacTransport);
      state.cov = BoundSymMatrix(
          std::get<FreeToBoundMatrix>(state.jacobian) *
          std::get<FreeSymMatrix>(state.cov) *
          std::get<FreeToBoundMatrix>(state.jacobian).transpose());
    } else {
      state.jacobian = FreeMatrix(state.jacTransport);
      state.cov =
          FreeSymMatrix(std::get<FreeMatrix>(state.jacobian) *
                        std::get<FreeSymMatrix>(state.cov) *
                        std::get<FreeMatrix>(state.jacobian).transpose());
    }
  }

  state.jacTransport = FreeMatrix::Identity();
  state.derivative = FreeVector::Zero();
  if (toLocal)
    reinitializeJacToGlobal(state);
  else
    state.jacToGlobal = std::nullopt;
}

}  // namespace detail
}  // namespace Acts