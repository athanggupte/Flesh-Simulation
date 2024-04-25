cmake_minimum_required(VERSION 3.19)
project(FleshSimulation LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# OpenMP
find_package(OpenMP)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)
igl_include(imgui)
igl_include(copyleft tetgen)

# Library with core functionality
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
add_library(flesh ${SOURCES} ${HEADERS})
target_include_directories(flesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(flesh
PUBLIC 
  igl::core
  igl::opengl
  igl::glfw

PRIVATE
  OpenMP::OpenMP_CXX
  )


# Driver project
add_executable(driver
  main.cpp
  tetrahedralize.cpp
  getopt.c)
target_include_directories(driver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link igl (and the glfw module) to your project
target_link_libraries(driver PUBLIC 
  igl::core
  igl::opengl
  igl::glfw
  igl::imgui
  igl_copyleft::tetgen
  flesh
  )
