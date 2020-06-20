// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <stdexcept>
#include <string>

#include "Acts/EventData/Measurement.hpp"

namespace FW {

/// Source link class for Alpide pixel hit.
///
/// The source link stores the measuremts, surface
///
class PixelSourceLink {
 public:
  PixelSourceLink(const Acts::Surface& surface, Acts::Vector2D values,
                  Acts::ActsSymMatrixD<2> cov)
      : m_values(values), m_cov(cov), m_surface(&surface) {}
  /// Must be default_constructible to satisfy SourceLinkConcept.
  PixelSourceLink() = default;
  PixelSourceLink(PixelSourceLink&&) = default;
  PixelSourceLink(const PixelSourceLink&) = default;
  PixelSourceLink& operator=(PixelSourceLink&&) = default;
  PixelSourceLink& operator=(const PixelSourceLink&) = default;

  constexpr const Acts::Surface& referenceSurface() const { return *m_surface; }

  Acts::FittableMeasurement<PixelSourceLink> operator*() const {
    return Acts::Measurement<PixelSourceLink, Acts::ParDef::eLOC_0,
                             Acts::ParDef::eLOC_1>{
        m_surface->getSharedPtr(), *this, m_cov, m_values[0], m_values[1]};
  }

 private:
  Acts::Vector2D m_values;
  Acts::ActsSymMatrixD<2> m_cov;
  // need to store pointers to make the object copyable
  const Acts::Surface* m_surface;
  friend bool operator==(const PixelSourceLink& lhs,
                         const PixelSourceLink& rhs) {
    return lhs.m_values.isApprox(rhs.m_values) and
           lhs.m_cov.isApprox(rhs.m_cov);
  }
};

}  // end of namespace FW
