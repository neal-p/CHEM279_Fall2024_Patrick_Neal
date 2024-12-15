file(GLOB SOURCES "*/*.cpp")
add_library(combined STATIC ${SOURCES})

target_include_directories(combined PUBLIC ${CMAKE_SOURCE_DIR}/include/)
target_link_libraries(combined PUBLIC ${LINK_LIBS})
