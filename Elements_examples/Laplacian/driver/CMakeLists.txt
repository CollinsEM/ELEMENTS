# Laplace target
set(Laplace_SRC_H ../include/header.h ../include/state.h ../include/utilities.h) # ../elements/elements.h ../elements/matar.h)
set(Laplace_SRC_CXX Pseudo_Laplacian.cpp)

add_executable(Laplace ${Laplace_SRC_CXX} ${Laplace_SRC_H})

target_link_libraries (Laplace elements)
target_link_libraries (Laplace input)

install(TARGETS Laplace DESTINATION bin)
# 

# Laplace_Kokkos target
if (KOKKOS)
    set(Laplace_Kokkos_SRC_H ../include/header.h ../include/state.h ../include/utilities.h) # ../elements/elements.h ../elements/matar.h)
    set(Laplace_Kokkos_SRC_CXX Pseudo_Laplacian_Kokkos.cpp)

    add_executable(Laplace_Kokkos ${Laplace_Kokkos_SRC_CXX} ${Laplace_Kokkos_SRC_H})

    set_target_properties(Laplace_Kokkos PROPERTIES COMPILE_DEFINITIONS HAVE_KOKKOS)

    if (CUDA)
        set_target_properties(Laplace_Kokkos PROPERTIES COMPILE_DEFINITIONS HAVE_CUDA)
    elseif (OPENMP)
        set_target_properties(Laplace_Kokkos PROPERTIES COMPILE_DEFINITIONS HAVE_OPENMP)
    endif()

    target_link_libraries (Laplace_Kokkos Kokkos::kokkos)
    target_link_libraries (Laplace_Kokkos kelements kinput)

    install(TARGETS Laplace_Kokkos DESTINATION bin)
endif()
#
