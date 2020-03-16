// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/ParameterTypes.hpp"

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
  // The following alias without prefix are
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
  // Distance-of-closest approach on a virtual perigee surface
  eLOC_D0 = eLOC_0,
  eLOC_Z0 = eLOC_1,
  ePHI = eBoundPhi,
  eTHETA = eBoundTheta,
  eQOP = eBoundSignedInverseP,
  eT = eBoundTime,
  BoundParsDim = eBoundParametersSize,
};

/// Underlying fundamental scalar type for bound track parameters.
using BoundParametersScalar = double;

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
  // For backward compatibility
  FreeParsDim = eFreeParametersSize,
};

/// Underlying fundamental scalar type for free track parameters.
using FreeParametersScalar = double;

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
  // for backward compatibility
  SpacePointDim = eSpacePointSize,
};

/// Underlying fundamental scalar type for space points.
using SpacePointScalar = double;

}  // namespace Acts
