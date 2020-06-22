// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TelescopeDetectorElement.h, Acts project, Generic Detector plugin
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/LineSurfaceStub.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;
class PlanarBounds;
class DiscBounds;
class ISurfaceMaterial;
class LineBounds;
}  // namespace Acts

namespace FW {

namespace Telescope {

/// @class TelescopeDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class TelescopeDetectorElement : public Acts::DetectorElementBase {
 public:
  /// Broadcast the ContextType
  using ContextType = Acts::GeometryContext;

  TelescopeDetectorElement() : Acts::DetectorElementBase() {}

  TelescopeDetectorElement(std::shared_ptr<const Acts::Transform3D> transform)
      : Acts::DetectorElementBase(), m_elementTransform(std::move(transform)) {}

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  TelescopeDetectorElement(
      std::shared_ptr<const Acts::Transform3D> transform,
      std::shared_ptr<const Acts::PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr)
      : Acts::DetectorElementBase(),
        m_elementTransform(std::move(transform)),
        m_elementThickness(thickness) {
    auto mutableSurface =
        Acts::Surface::makeShared<Acts::PlaneSurface>(pBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the line bounds for the line like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  TelescopeDetectorElement(
      std::shared_ptr<const Acts::Transform3D> transform,
      std::shared_ptr<const Acts::LineBounds> lBounds, double thickness,
      std::shared_ptr<const Acts::ISurfaceMaterial> material = nullptr)
      : Acts::DetectorElementBase(),
        m_elementTransform(std::move(transform)),
        m_elementThickness(thickness) {
    auto mutableSurface =
        Acts::Surface::makeShared<Acts::LineSurfaceStub>(lBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  ///  Destructor
  ~TelescopeDetectorElement() override { /*nop */
  }

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Acts::Transform3D& transform(
      const Acts::GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Acts::Surface& surface() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

 private:
  /// the transform for positioning in 3D space
  std::shared_ptr<const Acts::Transform3D> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Acts::Surface> m_elementSurface{nullptr};
  /// the element thickness
  double m_elementThickness{0.};
};

inline const Acts::Transform3D& TelescopeDetectorElement::transform(
    const Acts::GeometryContext& /*gctx*/) const {
  return *m_elementTransform;
}

inline const Acts::Surface& TelescopeDetectorElement::surface() const {
  return *m_elementSurface;
}

inline double TelescopeDetectorElement::thickness() const {
  return m_elementThickness;
}
}  // namespace Telescope
}  // namespace FW
