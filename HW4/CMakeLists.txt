add_executable(hw4 problem.cpp)
target_link_libraries(hw4 PRIVATE combined ${LINK_LIBS})
set_target_properties(hw4 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW4)

add_executable(scan scan_bond_distance.cpp)
target_link_libraries(scan PRIVATE combined ${LINK_LIBS})
set_target_properties(scan PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW4)


