
add_subdirectory(cycle_generation)

set(CURRENT_INCLUDE_DIR ${INCLUDE_DIR})
add_executable(main
        main.cpp
        Engine.cpp ${CURRENT_INCLUDE_DIR}/Engine.hpp
        ${CURRENT_INCLUDE_DIR}/EngineModeExecutor.hpp
        MatrixSizeExecutor.cpp ${CURRENT_INCLUDE_DIR}/MatrixSizeExecutor.hpp
        UpdateCountExecutor.cpp ${CURRENT_INCLUDE_DIR}/UpdateCountExecutor.hpp versioning.h.in ../headers/cycle_generation/CApi.h)

target_include_directories(main PUBLIC ${CURRENT_INCLUDE_DIR})

set_target_properties(main PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_BIN_DIRECTORY})

target_link_libraries(main PUBLIC cycle_generation OpenMP::OpenMP_CXX)


add_executable(c_api_test c_api_test.c)

target_include_directories(c_api_test PUBLIC ${CURRENT_INCLUDE_DIR})

set_target_properties(c_api_test PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${GLOBAL_BIN_DIRECTORY})

target_link_libraries(c_api_test PUBLIC cycle_generation)