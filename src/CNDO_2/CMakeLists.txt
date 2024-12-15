file(GLOB SOURCES "*.cpp")
add_library(CNDO_2 STATIC ${SOURCES})

target_include_directories(CNDO_2 PUBLIC ${CMAKE_SOURCE_DIR}/include/CNDO_2)
target_link_libraries(CNDO_2 PUBLIC ${LINK_LIBS})
