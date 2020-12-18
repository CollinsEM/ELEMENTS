cmake_minimum_required (VERSION 3.2)
project (Laplace)

set (CMAKE_CXX_STANDARD 17)

if (KOKKOS)
    #add_definitions(-DHAVE_KOKKOS=1)
  if (CUDA)
      #add_definitions(-DHAVE_CUDA=1)
    set(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
    set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE ON CACHE BOOL "")
  elseif (OPENMP)
      #add_definitions(-DHAVE_OPENMP=1)
    set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
  endif ()
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "")
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# --- custom targets: ---
INCLUDE( cmake/Modules/TargetDistclean.cmake OPTIONAL)

include_directories(./)
if (KOKKOS)
    include_directories(Kokkos)
endif()
include_directories(elements)
include_directories(include)


if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
 # set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

if (KOKKOS)
    add_subdirectory(Kokkos/kokkos)
    get_property(KOKKOS_INCLUDE_DIRS DIRECTORY Kokkos/kokkos/core/src PROPERTY INCLUDE_DIRECTORIES)
    include_directories(${KOKKOS_INCLUDE_DIRS})
endif()

add_subdirectory(elements)
add_subdirectory(src)
add_subdirectory(driver)
