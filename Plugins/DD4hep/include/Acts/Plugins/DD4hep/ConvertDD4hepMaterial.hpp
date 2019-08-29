// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <DD4hep/DetElement.h>
#include <DD4hep/DetFactoryHelper.h>
#include "Acts/Plugins/DD4hep/ActsExtension.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"

namespace Acts {

class ISurfaceMaterial;
class Layer;

/// Helper method to translate DD4hep material to Acts::ISurfaceMaterial
///
/// Thisis used to assign proto material to Cylinder Layers
///
/// @param detElement the DD4hep detector element for which this maerial is
/// assigned
/// @param loggingLevel is the output level for the conversion
///
/// @return a map of the identification string and a surface material
void addCylinderProtoMaterial(
    dd4hep::DetElement detElement, Layer& cylinderLayer,
    Logging::Level loggingLevel = Logging::Level::INFO);

/// Helper method to translate DD4hep material to Acts::ISurfaceMaterial
///
/// Thisis used to assign proto material to Disc Layers
///
/// @param detElement the DD4hep detector element for which this maerial is
/// assigned
/// @param loggingLevel is the output level for the conversion
///
/// @return a map of the identification string and a surface material
void addDiscProtoMaterial(dd4hep::DetElement detElement, Layer& discLayer,
                          Logging::Level loggingLevel = Logging::Level::INFO);

/// Helper method to be called for Cylinder and Disc Proto material
///
/// For both, cylinder and disc, the closed binning value is "binPhi"
///
/// @param aExtension the ActsExtension for the binning parameters
/// @param layer the Layer to assign the proto material
/// @param openBinning the open part of the binning (binR, binZ)
/// @param openBinValue the corresponding binning value
void addProtoMaterial(const ActsExtension& actsExtension, Layer& layer,
                      const std::string& openBinning,
                      Acts::BinningValue openBinVal);

/// Helper method that decorates an ActsExtension with proto material
/// description,
/// - it assigns bins for inner / representing / outer
///
/// @param x_layer the cylinder layer
/// @param actsExtension the extension that is augmented
/// @param materialOptions the placement options (inner / representing / outer)
/// @param binOptions the material binning options
void xml2LayerProtoMaterial(
    const xml_comp_t& x_layer, ActsExtension& actsExtension,
    const std::vector<std::string>& materialOptions,
    const std::pair<std::string, std::string>& binOptions);

/// Helper method that decorates an ActsExtension with proto material
/// description,
/// - it assigns bins for inner / representing / outer
///
/// @param x_layer the cylinder layer
/// @param actsExtension the extension that is augmented
void xml2CylinderProtoMaterial(const xml_comp_t& x_layer,
                               Acts::ActsExtension& actsExtension);

/// Helper method that decorates an ActsExtension with proto material
/// description,
/// - it assigns bins for inner / representing / outer
///
/// @param x_layer the disc layer
/// @param actsExtension the extension that is augmented
void xml2DiscProtoMaterial(const xml_comp_t& x_layer,
                           Acts::ActsExtension& actsExtension);

}  // namespace Acts