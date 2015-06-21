cmake_minimum_required (VERSION 2.6) 
project(maxwellDist)

set(CMAKE_CXX_FLAGS "-fPIC -g -Wall -Wno-write-strings")

file(GLOB src_file1 Particle*.cxx)
file(GLOB inc_file1 Particle*.h)
add_library(Particle SHARED ${src_file1} ${inc_file1})
target_link_libraries(Particle Vector)

file(GLOB src_file2 Vector*.cxx)
file(GLOB inc_file2 Vector*.h)
add_library(Vector SHARED ${src_file2} ${inc_file2})
target_link_libraries(Vector)

file(GLOB src_file3 maxwellDist*.cxx)
file(GLOB inc_file3 maxwellDist*.h)
add_library(MaxwellDist SHARED ${src_file3} ${inc_file3})
target_link_libraries(MaxwellDist Particle Vector)

ADD_EXECUTABLE(main.exe "main.cxx")
TARGET_LINK_LIBRARIES(main.exe Vector Particle MaxwellDist)
