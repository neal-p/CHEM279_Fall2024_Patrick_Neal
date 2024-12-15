add_executable(hw3 problem.cpp)
target_link_libraries(hw3 PRIVATE combined ${LINK_LIBS})
set_target_properties(hw3 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW3)

add_executable(bond_energy bond_energy.cpp)
target_link_libraries(bond_energy PRIVATE combined ${LINK_LIBS})
set_target_properties(bond_energy PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW3)

add_executable(basic_rxn basic_rxn.cpp)
target_link_libraries(basic_rxn PRIVATE combined ${LINK_LIBS})
set_target_properties(basic_rxn PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW3)

