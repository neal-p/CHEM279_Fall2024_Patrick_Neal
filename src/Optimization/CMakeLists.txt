file(GLOB SOURCES "*.cpp")
add_library(Optimization STATIC ${SOURCES})

target_include_directories(Optimization PUBLIC ${CMAKE_SOURCE_DIR}/include/Optimization)
target_link_libraries(Optimization PUBLIC ${LINK_LIBS})
