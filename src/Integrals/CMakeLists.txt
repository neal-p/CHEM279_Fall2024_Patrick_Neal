file(GLOB SOURCES "*.cpp")
add_library(Integrals STATIC ${SOURCES})

target_include_directories(Integrals PUBLIC ${CMAKE_SOURCE_DIR}/include/Integrals)
target_link_libraries(Integrals PUBLIC ${LINK_LIBS})
