cmake_minimum_required(VERSION 3.0)

project(drop)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(BasicConfig)

cond_find_package(mia REQUIRED)
cond_find_package(mrfopt REQUIRED)

find_package(TBB REQUIRED)

file(GLOB_RECURSE SRC "src/*.cpp")
file(GLOB_RECURSE HDR "include/*.h")
list(APPEND SRC ${HDR})

add_library(${PROJECT_NAME} STATIC ${SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_PDB_NAME_DEBUG "${PROJECT_NAME}${CMAKE_DEBUG_POSTFIX}")

target_include_directories(${PROJECT_NAME} PUBLIC "include")
target_include_directories(${PROJECT_NAME} PUBLIC ${TBB_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC "mia")
target_link_libraries(${PROJECT_NAME} PUBLIC "mrfopt")
target_link_libraries(${PROJECT_NAME} PUBLIC ${TBB_LIBRARIES})

export(PACKAGE drop)
export(TARGETS drop
  FILE "${PROJECT_BINARY_DIR}/dropConfig.cmake"
)
