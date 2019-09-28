cmake_minimum_required(VERSION 3.0)

project(itkio)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(BasicConfig)

cond_find_package(mia REQUIRED)

find_package(ITK REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)

include(${ITK_USE_FILE})

file(GLOB_RECURSE SRC "src/*.cpp")
file(GLOB_RECURSE HDR "include/*.h")
list(APPEND SRC ${HDR})

add_library(${PROJECT_NAME} STATIC ${SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_PDB_NAME_DEBUG "${PROJECT_NAME}${CMAKE_DEBUG_POSTFIX}")

target_include_directories(${PROJECT_NAME} PUBLIC "include")
target_include_directories(${PROJECT_NAME} PUBLIC ${ITK_INCLUDE_DIRS})
target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} PUBLIC "mia")
target_link_libraries(${PROJECT_NAME} PUBLIC ${ITK_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_SYSTEM_LIBRARY})
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_FILESYSTEM_LIBRARY})

export(PACKAGE itkio)
export(TARGETS itkio
  FILE "${PROJECT_BINARY_DIR}/itkioConfig.cmake"
)
