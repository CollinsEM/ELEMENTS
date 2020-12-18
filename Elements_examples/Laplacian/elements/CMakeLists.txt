add_library(elements STATIC elements.cpp elements.h matar.h)

# if we are using Kokkos
if (KOKKOS)
  add_library(kelements elements.cpp elements.h matar.h)
  set_target_properties(kelements PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)

  if (CUDA)
      set_target_properties(kelements PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
  elseif (OPENMP)
      set_target_properties(kelements PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
  endif()

  target_link_libraries(kelements Kokkos::kokkos)
endif()
