add_executable(hw1_1 problem1.cpp)
target_link_libraries(hw1_1 PRIVATE Simulation Energy ${LINK_LIBS})
set_target_properties(hw1_1 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW1)

add_executable(hw1_2 problem2.cpp)
target_link_libraries(hw1_2 PRIVATE Simulation Energy Gradients ${LINK_LIBS})
set_target_properties(hw1_2 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW1)

add_executable(hw1_3 problem3.cpp)
target_link_libraries(hw1_3 PRIVATE Simulation Energy Gradients Optimization ${LINK_LIBS})
set_target_properties(hw1_3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW1)

