include_directories(Kokkos)
add_subdirectory(Kokkos/kokkos)
get_property(KOKKOS_INCLUDE_DIRS DIRECTORY Kokkos/kokkos PROPERTY INCLUDE_DIRECTORIES)
include_directories(${KOKKOS_INCLUDE_DIRS})

add_library(elements elements.cpp elements.h matar.h)
