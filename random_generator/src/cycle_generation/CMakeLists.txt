
add_subdirectory(update_generation)

set(CURRENT_INCLUDE_DIR ${INCLUDE_DIR}/cycle_generation)
add_library(cycle_generation
        ${CURRENT_INCLUDE_DIR}/Cycle.hpp
        Matrix.cpp ${CURRENT_INCLUDE_DIR}/Matrix.hpp
        hdf5CycleOutputStream.cpp ${CURRENT_INCLUDE_DIR}/hdf5CycleOutputStream.hpp
        SizedBasedCycleStream.cpp ${CURRENT_INCLUDE_DIR}/SizedBasedCycleStream.hpp
        CycleGenerator.cpp ${CURRENT_INCLUDE_DIR}/CycleGenerator.hpp
        CApi.cpp)
target_include_directories(cycle_generation PUBLIC ${CURRENT_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS})
target_link_libraries(cycle_generation PUBLIC update_generation ${HDF5_LIBRARIES}
        m ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} lapacke)