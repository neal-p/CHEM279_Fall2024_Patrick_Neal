file(GLOB SOURCES "*.cpp")
add_library(Energy STATIC ${SOURCES})

target_include_directories(Energy PUBLIC ${CMAKE_SOURCE_DIR}/include/Energy)
target_link_libraries(Energy PUBLIC ${LINK_LIBS})
