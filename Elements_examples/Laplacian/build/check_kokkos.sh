./driver/Laplace_Kokkos ../meshes/mesh_2x2x2.geo >& ../ktest/lap.kokkos
./driver/Laplace ../meshes/mesh_2x2x2.geo >& ../ktest/lap.silver
vimdiff ../ktest/lap.kokkos ../ktest/lap.silver
