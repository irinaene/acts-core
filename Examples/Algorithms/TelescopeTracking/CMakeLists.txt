add_library(
  ActsExamplesTelescopeTracking SHARED
  src/TelescopeTrackingAlgorithm.cpp
  src/TelescopeTrackingAlgorithmFitterFunction.cpp)
target_include_directories(
  ActsExamplesTelescopeTracking
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesTelescopeTracking
  PUBLIC
    ActsCore
    ActsExamplesFramework ActsExamplesMagneticField
    Boost::program_options)

install(
  TARGETS ActsExamplesTelescopeTracking
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
