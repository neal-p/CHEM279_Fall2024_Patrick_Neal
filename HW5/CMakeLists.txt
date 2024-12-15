  add_executable(hw5 problem.cpp)
  target_link_libraries(hw5 PRIVATE combined ${LINK_LIBS})
  set_target_properties(hw5 PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/HW5)
