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
        # include 下的所有头文件 (仅为了视图显示)
)
# 静态链接标准库，避免依赖 DLL
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")