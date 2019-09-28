cmake_minimum_required(VERSION 2.8.12)

project(mrfopt)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(BasicConfig)

file(GLOB_RECURSE MRFOPT_SRC "src/*.cpp")
file(GLOB_RECURSE MRFOPT_HDR "include/*.h")
file(GLOB_RECURSE MRFOPT_HPP "include/*.hpp")
file(GLOB_RECURSE MRFOPT_EXT_SRC "ext/*.cpp")
file(GLOB_RECURSE MRFOPT_EXT_HDR "ext/*.h")
list(APPEND MRFOPT_SRC ${MRFOPT_HDR})
list(APPEND MRFOPT_SRC ${MRFOPT_HPP})
list(APPEND MRFOPT_SRC ${MRFOPT_EXT_SRC})
list(APPEND MRFOPT_SRC ${MRFOPT_EXT_HDR})

add_library(${PROJECT_NAME} STATIC ${MRFOPT_SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_PDB_NAME_DEBUG "${PROJECT_NAME}${CMAKE_DEBUG_POSTFIX}")

target_include_directories(${PROJECT_NAME} PUBLIC "include")
target_include_directories(${PROJECT_NAME} PUBLIC "ext/elc")
target_include_directories(${PROJECT_NAME} PUBLIC "ext/qpbo")

export(PACKAGE mrfopt)
export(TARGETS mrfopt
  FILE "${PROJECT_BINARY_DIR}/mrfoptConfig.cmake"
)
