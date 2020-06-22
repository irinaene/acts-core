add_executable(ActsRecTelescopeTracks ./RecTelescopeTracks.cpp)
target_link_libraries(ActsRecTelescopeTracks  
  PRIVATE
  ActsExamplesFramework
  ActsExamplesCommon
  ActsExamplesTelescopeTracking
  ActsExamplesDetectorTelescope
  ActsExamplesMagneticField
  ActsExamplesIoCsv
  ActsExamplesIoObj
  ActsExamplesIoPerformance
  altel_acts
  )


install(TARGETS ActsRecTelescopeTracks
  # EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin      COMPONENT runtime
  LIBRARY       DESTINATION lib64    COMPONENT runtime
  ARCHIVE       DESTINATION lib64    COMPONENT devel
  PUBLIC_HEADER DESTINATION include  COMPONENT devel
  RESOURCE      DESTINATION resource COMPONENT runtime
  )
