file(GLOB SOURCES "*.cpp")
add_library(Gaussian STATIC ${SOURCES})

target_include_directories(Gaussian PUBLIC ${CMAKE_SOURCE_DIR}/include/Gaussian)
target_link_libraries(Gaussian PUBLIC ${LINK_LIBS})
