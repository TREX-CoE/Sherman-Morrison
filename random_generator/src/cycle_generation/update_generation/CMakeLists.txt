
set(CURRENT_INCLUDE_DIR ${INCLUDE_DIR}/cycle_generation/update_generation)
add_library(update_generation
        UpdateDescriptor.cpp ${CURRENT_INCLUDE_DIR}/UpdateDescriptor.hpp
        UpdateMatrix.cpp ${CURRENT_INCLUDE_DIR}/UpdateMatrix.hpp
        UpdatePlan.cpp ${CURRENT_INCLUDE_DIR}/UpdatePlan.hpp
        ${CURRENT_INCLUDE_DIR}/UpdatePlanner.hpp
        UpdateMatrixGenerator.cpp ${CURRENT_INCLUDE_DIR}/UpdateGenerator.hpp
        UpdateGenerator.cpp ${CURRENT_INCLUDE_DIR}/UpdateGenerator.hpp
        ChainedUpdatePlanner.cpp ${CURRENT_INCLUDE_DIR}/ChainedUpdatePlanner.hpp
        )
set_target_properties(update_generation PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY ${GLOBAL_LIB_DIRECTORY} LIBRARY_OUTPUT_DIRECTORY ${GLOBAL_LIB_DIRECTORY}
        )

target_include_directories(update_generation PUBLIC ${CURRENT_INCLUDE_DIR} ${INCLUDE_DIR})