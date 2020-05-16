cmake_minimum_required (VERSION 3.1)

project (RGBDPlaneDetection)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fpermissive")
set(CMAKE_CXX_COMPILER /usr/bin/g++)
set(CMAKE_BUILD_TYPE "Release")

# put FindGLFW.cmake in the source directory of CMakeLists.txt
find_package(OpenCV 3.4 REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS} /usr/local/include/eigen3)
include_directories(include/peac)
include_directories(include/MRF2.2)

link_directories(include/MRF2.2)

file(GLOB sources "*.cpp")

add_executable(RGBDPlaneDetection ${sources})

target_link_libraries(RGBDPlaneDetection ${OpenCV_LIBS} libMRF.a)