file(GLOB SOURCES "*.cpp")
add_library(Simulation STATIC ${SOURCES})

target_include_directories(Simulation PUBLIC ${CMAKE_SOURCE_DIR}/include/Simulation)
target_link_libraries(Simulation PUBLIC ${LINK_LIBS})
