// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <type_traits>

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/ParameterTypes.hpp"

// The user can override the (track) parameter ordering and underlying scalar
// type. If the variable is defined, it must point to a header file that
// contains the same enum and type definitions for bound and free track
// parameters as well as space points as given below.
#ifdef ACTS_PARAMETER_DEFINITIONS_HEADER
#include ACTS_PARAMETER_DEFINITIONS_HEADER
#else
namespace Acts {

/// Components of a bound track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum BoundParametersIndices : unsigned int {
  eBoundLoc0 = 0,
  eBoundLoc1 = 1,
  eBoundPhi = 2,
  eBoundTheta = 3,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  eBoundSignedInverseP = 4,
  eBoundTime = 5,
  // Last uninitialized value contains the total number of components
  eBoundParametersSize,
  // The following aliases without prefix exist for historical reasons
  // Generic spatial coordinates on the local surface
  eLOC_0 = eBoundLoc0,
  eLOC_1 = eBoundLoc1,
  // Spatial coordinates on a disk in polar coordinates
  eLOC_R = eLOC_0,
  eLOC_PHI = eLOC_1,
  // Spatial coordinates on a disk in Cartesian coordinates
  eLOC_X = eLOC_0,
  eLOC_Y = eLOC_1,
  // Spatial coordinates on a cylinder
  eLOC_RPHI = eLOC_0,
  eLOC_Z = eLOC_1,
  // Closest approach coordinates on a virtual perigee surface
  eLOC_D0 = eLOC_0,
  eLOC_Z0 = eLOC_1,
  ePHI = eBoundPhi,
  eTHETA = eBoundTheta,
  eQOP = eBoundSignedInverseP,
  eT = eBoundTime,
};

/// Components of a free track parameters vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum FreeParametersIndices : unsigned int {
  // Cartesian position
  // The spatial position components must be stored as one continous block.
  eFreePos0 = 0u,
  eFreePos1 = eFreePos0 + 1u,
  eFreePos2 = eFreePos0 + 2u,
  // Time
  eFreeTime = 3u,
  // Cartesian (unit) direction
  // The direction components must be stored as one continous block.
  eFreeDir0 = 4u,
  eFreeDir1 = eFreeDir0 + 1u,
  eFreeDir2 = eFreeDir0 + 2u,
  // Global inverse-momentum-like parameter, i.e. q/p or 1/p
  eFreeSignedInverseP = 7u,
  // Last uninitialized value contains the total number of components
  eFreeParametersSize,
};

/// Components of a space point vector.
///
/// To be used to access components by named indices instead of just numbers.
/// This must be a regular `enum` and not a scoped `enum class` to allow
/// implicit conversion to an integer. The enum value are thus visible directly
/// in `namespace Acts` and are prefixed to avoid naming collisions.
enum SpacePointIndices : unsigned int {
  // For spatial points
  // The spatial position components must be stored as one continous block.
  eSpX = 0u,
  eSpY = eSpX + 1u,
  eSpZ = eSpX + 2u,
  eSpT = 3u,
  // Alias to allow clear code when accessing momentum vectors
  eSpPx = eSpX,
  eSpPy = eSpY,
  eSpPz = eSpZ,
  eSpE = eSpT,
  // Last uninitialized value contains the total number of components
  eSpacePointSize,
};

/// Underlying fundamental scalar type for bound track parameters.
using BoundParametersScalar = double;
/// Underlying fundamental scalar type for free track parameters.
using FreeParametersScalar = double;
/// Underlying fundamental scalar type for space points.
using SpacePointScalar = double;

}  // namespace Acts
#endif

namespace Acts {

// Ensure bound track parameters definition is valid.
static_assert(std::is_enum_v<BoundParametersIndices>,
              "'BoundParametersIndices' must be an enum type");
static_assert(std::is_convertible_v<BoundParametersIndices, size_t>,
              "'BoundParametersIndices' must be convertible to size_t");
static_assert(2 <= BoundParametersIndices::eBoundParametersSize,
              "Bound track parameters must have at least two components");
static_assert(std::is_floating_point_v<BoundParametersScalar>,
              "'BoundParametersScalar' must be a floating point type");

// Ensure free track parameters definition is valid.
static_assert(std::is_enum_v<FreeParametersIndices>,
              "'FreeParametersIndices' must be an enum type");
static_assert(std::is_convertible_v<FreeParametersIndices, size_t>,
              "'FreeParametersIndices' must be convertible to size_t");
static_assert(6 <= FreeParametersIndices::eFreeParametersSize,
              "Free track parameters must have at least six components");
static_assert(std::is_floating_point_v<FreeParametersScalar>,
              "'FreeParametersScalar' must be a floating point type");

// Ensure space point definition is valid.
static_assert(std::is_enum_v<SpacePointIndices>,
              "'SpacePointIndices' is not an enum type");
static_assert(std::is_convertible_v<SpacePointIndices, size_t>,
              "'SpacePointIndices' is not convertible to size_t");
static_assert(3 <= SpacePointIndices::eSpacePointSize,
              "Space points must have at least three components");
static_assert(std::is_floating_point_v<SpacePointScalar>,
              "'SpacePointScalar' must be a floating point type");

// Ensure bound track parameter components/ indices are consistently defined.
static_assert(eLOC_0 != eLOC_1, "Local parameters must be differents");
static_assert(eLOC_R == eLOC_0 or eLOC_R == eLOC_1,
              "Local radius must be a local parameter");
static_assert(eLOC_PHI == eLOC_0 or eLOC_PHI == eLOC_1,
              "Local phi must be a local parameter");
static_assert(eLOC_RPHI == eLOC_0 or eLOC_RPHI == eLOC_1,
              "Local r*phi must be a local parameter");
static_assert(eLOC_Z == eLOC_0 or eLOC_Z == eLOC_1,
              "Local z must be a local parameter");
static_assert(eLOC_X == eLOC_0 or eLOC_X == eLOC_1,
              "Local x must be a local parameter");
static_assert(eLOC_Y == eLOC_0 or eLOC_Y == eLOC_1,
              "Local y must be a local parameter");
static_assert(eLOC_D0 == eLOC_0 or eLOC_D0 == eLOC_1,
              "D0 must be a local parameter");
static_assert(eLOC_Z0 == eLOC_0 or eLOC_Z0 == eLOC_1,
              "Z0 must be a local parameter");

// Ensure free track parameter components/ indices are consistently defined.
static_assert(eFreePos1 == eFreePos0 + 1u, "Position must be continous");
static_assert(eFreePos2 == eFreePos0 + 2u, "Position must be continous");
static_assert(eFreeDir1 == eFreeDir0 + 1u, "Direction must be continous");
static_assert(eFreeDir2 == eFreeDir0 + 2u, "Direction must be continous");

// Ensure space point components/ indices are consistently defined.
static_assert(eSpY == eSpX + 1u, "Position must be continous");
static_assert(eSpZ == eSpX + 2u, "Position must be continous");
static_assert(eSpX == eSpPx, "Position and momentum must be consistent");
static_assert(eSpY == eSpPy, "Position and momentum must be consistent");
static_assert(eSpZ == eSpPz, "Position and momentum must be consistent");
static_assert(eSpT == eSpE, "Position and momentum must be consistent");

/// Bound track parameters traits for component constrains.
template <BoundParametersIndices>
struct BoundParametersTraits;
template <>
struct BoundParametersTraits<BoundParametersIndices::eLOC_0> {
  using type = local_parameter;
};
template <>
struct BoundParametersTraits<BoundParametersIndices::eLOC_1> {
  using type = local_parameter;
};
template <>
struct BoundParametersTraits<BoundParametersIndices::ePHI> {
  static constexpr double pMin() { return -M_PI; }
  static constexpr double pMax() { return M_PI; }
  using type = cyclic_parameter<double, pMin, pMax>;
};
template <>
struct BoundParametersTraits<BoundParametersIndices::eTHETA> {
  static constexpr double pMin() { return 0; }
  static constexpr double pMax() { return M_PI; }
  using type = bound_parameter<double, pMin, pMax>;
};
template <>
struct BoundParametersTraits<BoundParametersIndices::eQOP> {
  using type = unbound_parameter;
};
template <>
struct BoundParametersTraits<BoundParametersIndices::eT> {
  using type = unbound_parameter;
};

// The following matrix and vector types are automatically derived from the
// indices enums and scalar typedefs.

// Matrix and vector types related to bound track parameters.

using BoundVector = ActsVector<BoundParametersScalar, eBoundParametersSize>;
using BoundRowVector =
    ActsRowVector<BoundParametersScalar, eBoundParametersSize>;
using BoundMatrix = ActsMatrix<BoundParametersScalar, eBoundParametersSize,
                               eBoundParametersSize>;
using BoundSymMatrix =
    ActsSymMatrix<BoundParametersScalar, eBoundParametersSize>;

// Matrix and vector types related to free track parameters.

using FreeVector = ActsVector<FreeParametersScalar, eFreeParametersSize>;
using FreeRowVector = ActsRowVector<FreeParametersScalar, eFreeParametersSize>;
using FreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eFreeParametersSize>;
using FreeSymMatrix = ActsSymMatrix<FreeParametersScalar, eFreeParametersSize>;

// Matrix and vector types related to space points.

using SpacePointVector = ActsVector<SpacePointScalar, eSpacePointSize>;
using SpacePointRowVector = ActsRowVector<SpacePointScalar, eSpacePointSize>;
using SpacePointSymMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eSpacePointSize>;
using SpacePointSymMatrix = ActsSymMatrix<SpacePointScalar, eSpacePointSize>;

// Mapping to bound track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// bound track parameters. Thus, the bound parameters scalar type is sufficient
// to retain accuracy.

using FreeToBoundMatrix = ActsMatrix<BoundParametersScalar,
                                     eBoundParametersSize, eFreeParametersSize>;
using SpacePointToBoundMatrix =
    ActsMatrix<BoundParametersScalar, eBoundParametersSize, eSpacePointSize>;

// Mapping to free track parameters.
//
// Assumes that matrices represent maps from another space into the space of
// free track parameters. Thus, the free parameters scalar type is sufficient
// to retain accuracy.

using BoundToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eBoundParametersSize>;
using SpacePointToFreeMatrix =
    ActsMatrix<FreeParametersScalar, eFreeParametersSize, eSpacePointSize>;

// Mapping to space points.
//
// Assumes that matrices represent maps from another space into the space point
// space. Thus, the space point scalar type is sufficient to retain accuracy.

using BoundToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eBoundParametersSize>;
using FreeToSpacePointMatrix =
    ActsMatrix<SpacePointScalar, eSpacePointSize, eFreeParametersSize>;

// For backward compatibility. New code must use the more explicit
// `BoundParameters{Indices,Scalar,Traits}...` types.
using ParDef = BoundParametersIndices;
using ParID_t = BoundParametersIndices;
using ParValue_t = BoundParametersScalar;
template <BoundParametersIndices idx>
using par_type = BoundParametersTraits<idx>;
template <BoundParametersIndices idx>
using par_type_t = typename BoundParametersTraits<idx>::type;

}  // namespace Acts
