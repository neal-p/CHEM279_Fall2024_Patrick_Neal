file(GLOB SOURCES "*.cpp")
add_library(Math STATIC ${SOURCES})

target_include_directories(Math PUBLIC ${CMAKE_SOURCE_DIR}/include/Math)
target_link_libraries(Math PUBLIC ${LINK_LIBS})
