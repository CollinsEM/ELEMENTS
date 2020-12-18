add_library(input input.cpp)

# if we are using Kokkos
if (KOKKOS)
  add_library(kinput input.cpp)
  set_target_properties(kinput PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)

  if (CUDA)
      set_target_properties(kinput PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
  elseif (OPENMP)
      set_target_properties(kinput PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
  endif()

  target_link_libraries(input Kokkos::kokkos)
endif()
