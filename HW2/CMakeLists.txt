add_executable(hw2_1 problem1.cpp)
target_link_libraries(hw2_1 PRIVATE Simulation Integrands Integrals ${LINK_LIBS})
set_target_properties(hw2_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW2)

add_executable(hw2_2 problem2.cpp)
target_link_libraries(hw2_2 PRIVATE Simulation Integrands Integrals Math ${LINK_LIBS})
set_target_properties(hw2_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW2)

