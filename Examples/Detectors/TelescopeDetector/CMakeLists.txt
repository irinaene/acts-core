add_library(
  ActsExamplesDetectorTelescope SHARED
  src/TelescopeDetector.cpp)
target_include_directories(
  ActsExamplesDetectorTelescope
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesDetectorTelescope
  PUBLIC ActsCore ActsIdentificationPlugin ActsDigitizationPlugin)
target_link_libraries(
  ActsExamplesDetectorTelescope
  PUBLIC ActsExamplesFramework ActsExamplesDetectorsCommon)

install(
  TARGETS ActsExamplesDetectorTelescope
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
