// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <vector>

#include "ACTFW/TelescopeDetector/TelescopeDetectorElement.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/PassiveLayerBuilder.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {
class TrackingGeometry;
}

namespace FW {
namespace Telescope {

/// Global method to build the telescope tracking geometry
///
/// @tparam detector_element_t is the actual type of the detector
/// element, each derivative of a TelescopeDetectorElement can be used
///
/// @param gctx is the detector element dependent geometry context
/// @param detectorStore is the store for the detector element
/// @param matDecorator is an optional decorator for the material
template <typename detector_element_t>
std::unique_ptr<const Acts::TrackingGeometry> buildDetector(
    const typename detector_element_t::ContextType& gctx,
    std::vector<std::vector<std::shared_ptr<detector_element_t>>>&
    /*detectorStore*/,
    std::shared_ptr<const Acts::IMaterialDecorator> matDecorator = nullptr) {
  using namespace Acts::UnitLiterals;

  // Construct the rotation
  Acts::RotationMatrix3D rotation = Acts::RotationMatrix3D::Identity();
  double rotationAngle = 90_degree;
  Acts::Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
  Acts::Vector3D yPos(0., 1., 0.);
  Acts::Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
  rotation.col(0) = xPos;
  rotation.col(1) = yPos;
  rotation.col(2) = zPos;

  // Boundaries of the surfaces (ALPIDE SIZE: 30.81896 * 13.76256_mm)
  const auto rBounds = std::make_shared<const Acts::RectangleBounds>(
      Acts::RectangleBounds(30.81896_mm, 13.76256_mm));

  // Material of the surfaces
  Acts::MaterialProperties matProp(95.7, 465.2, 28.03, 14., 2.32e-3, 100_um);
  const auto surfaceMaterial =
      std::make_shared<Acts::HomogeneousSurfaceMaterial>(matProp);

  // Set translation vectors
  std::vector<Acts::Vector3D> translations;
  translations.reserve(6);
  translations.push_back({-95_mm, 0., 0.});
  translations.push_back({-57_mm, 0., 0.});
  translations.push_back({-19_mm, 0., 0.});
  translations.push_back({19_mm, 0., 0.});
  translations.push_back({57_mm, 0., 0.});
  translations.push_back({95_mm, 0., 0.});

  // Construct layer configs
  std::vector<Acts::CuboidVolumeBuilder::LayerConfig> lConfs;
  lConfs.reserve(6);
  unsigned int i;
  for (i = 0; i < translations.size(); i++) {
    Acts::CuboidVolumeBuilder::SurfaceConfig sConf;
    sConf.position = translations[i];
    sConf.rotation = rotation;
    sConf.rBounds = rBounds;
    sConf.surMat = surfaceMaterial;
    // The thickness to construct the associated detector element
    sConf.thickness = 1._um;
    sConf.detElementConstructor =
        [](std::shared_ptr<const Acts::Transform3D> trans,
           std::shared_ptr<const Acts::RectangleBounds> bounds,
           double thickness) {
          return new TelescopeDetectorElement(trans, bounds, thickness);
        };
    Acts::CuboidVolumeBuilder::LayerConfig lConf;
    lConf.surfaceCfg = sConf;
    lConfs.push_back(lConf);
  }

  // Construct volume config
  Acts::CuboidVolumeBuilder::VolumeConfig vConf;
  vConf.position = {0., 0., 0.};
  vConf.length = {1.2_m, 1._m, 1._m};
  vConf.layerCfg = lConfs;
  vConf.name = "Tracker";

  // Construct volume builder config
  Acts::CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {1.2_m, 1._m, 1._m};
  conf.volumeCfg = {vConf};  // one volume

  // Build detector
  Acts::CuboidVolumeBuilder cvb;
  cvb.setConfig(conf);
  Acts::TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto& vb) {
        return cvb.trackingVolume(context, inner, vb);
      });
  Acts::TrackingGeometryBuilder tgb(tgbCfg);
  auto detector = tgb.trackingGeometry(gctx);

  /// return the tracking geometry
  return detector;
}

}  // end of namespace Telescope
}  // end of namespace FW
