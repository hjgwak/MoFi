cmake_minimum_required(VERSION 3.7)
project(MoFi)

set(CMAKE_CXX_STANDARD 11)

add_library(Sequence Sequence.cpp)
add_library(SequenceReader SequenceReader.cpp Sequence.cpp)

set(SOURCE_FILES main.cpp)
add_executable(MoFi ${SOURCE_FILES})

target_link_libraries(MoFi Sequence SequenceReader)