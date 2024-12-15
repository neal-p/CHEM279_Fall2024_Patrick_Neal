cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANRARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


project(
            CHEM279
            VERSION 1.0
            LANGUAGES CXX
  )


# Global library includes
include_directories(${CMAKE_SOURCE_DIR}/include)

# Eigen includes
include_directories(${CMAKE_SOURCE_DIR}/external)

# Blas and Lapack
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
set(LINK_LIBS ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#####################################################################################

# Common libraries
# will be developed for re-use in multiple  hws
#add_subdirectory(src/Simulation)
#add_subdirectory(src/Energy)
#add_subdirectory(src/Gradients)
#add_subdirectory(src/Optimization)
#add_subdirectory(src/Gaussian)
#add_subdirectory(src/Integrals)
#add_subdirectory(src/Math)
add_subdirectory(src)

#####################################################################################


# Homework EXECUTABLES
#add_subdirectory(HW1)
#add_subdirectory(HW2)
#add_subdirectory(HW3)
#add_subdirectory(HW4)
add_subdirectory(HW5)

