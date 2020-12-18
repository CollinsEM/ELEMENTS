#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Kokkos::kokkoscore" for configuration "RelWithDebInfo"
set_property(TARGET Kokkos::kokkoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Kokkos::kokkoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib64/libkokkoscore.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Kokkos::kokkoscore )
list(APPEND _IMPORT_CHECK_FILES_FOR_Kokkos::kokkoscore "${_IMPORT_PREFIX}/lib64/libkokkoscore.a" )

# Import target "Kokkos::kokkoscontainers" for configuration "RelWithDebInfo"
set_property(TARGET Kokkos::kokkoscontainers APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(Kokkos::kokkoscontainers PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib64/libkokkoscontainers.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Kokkos::kokkoscontainers )
list(APPEND _IMPORT_CHECK_FILES_FOR_Kokkos::kokkoscontainers "${_IMPORT_PREFIX}/lib64/libkokkoscontainers.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
