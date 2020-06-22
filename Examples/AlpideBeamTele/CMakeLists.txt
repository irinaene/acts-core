add_subdirectory(rapidjson)

aux_source_directory(src LIB_SRC)
add_library(ActsExamplesAlpideBeamTele SHARED ${LIB_SRC})
target_include_directories(ActsExamplesAlpideBeamTele
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  )


target_link_libraries(
  ActsExamplesAlpideBeamTele
  
  PUBLIC
  Boost::program_options
  Threads::Threads

  ActsCore
  ActsDigitizationPlugin
  ActsIdentificationPlugin
  
  ActsExamplesFramework
  ActsExamplesDetectorsCommon
  ActsExamplesFramework
  ActsExamplesPropagation
  ActsExamplesTruthTracking
  ActsExamplesMagneticField
  
  ActsExamplesIoPerformance
  ActsExamplesIoRoot
  ActsExamplesIoObj
  ActsExamplesIoCsv
  
  PRIVATE
  ROOT::Core
  ROOT::Hist
  ROOT::Tree
  )

# set(LIB_PUBLIC_HEADERS include/xx.hh)
# set_target_properties(aidatlu PROPERTIES PUBLIC_HEADER "${LIB_PUBLIC_HEADERS}")

install(TARGETS ActsExamplesAlpideBeamTele
  # EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib64    COMPONENT runtime
  ARCHIVE       DESTINATION lib64    COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )

add_subdirectory(exe)





