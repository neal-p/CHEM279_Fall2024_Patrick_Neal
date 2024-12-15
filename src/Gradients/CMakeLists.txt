file(GLOB SOURCES "*.cpp")
add_library(Gradients STATIC ${SOURCES})

target_include_directories(Gradients PUBLIC ${CMAKE_SOURCE_DIR}/include/Gradients)
target_link_libraries(Gradients PUBLIC ${LINK_LIBS})
