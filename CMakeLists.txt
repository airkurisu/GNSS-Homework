cmake_minimum_required(VERSION 3.20)
project(mygnss10)
set(CMAKE_CXX_STANDARD 14)
include_directories(include)
file(GLOB SRC_FILES "src/*.cpp")
file(GLOB HEADER_FILES "include/*.h")
add_executable(mygnss10
        main.cpp
        ${SRC_FILES}
        ${HEADER_FILES}
)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
