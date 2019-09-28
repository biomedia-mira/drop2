cmake_minimum_required(VERSION 3.0)

project(mia_resample)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(BasicConfig)

cond_find_package(mia REQUIRED)
cond_find_package(itkio REQUIRED)

find_package(Boost COMPONENTS program_options REQUIRED)

file(GLOB_RECURSE SRC "src/*.cpp")
file(GLOB_RECURSE HDR "include/*.h")
list(APPEND SRC ${HDR})

add_executable(${PROJECT_NAME} ${SRC})
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_PDB_NAME_DEBUG "${PROJECT_NAME}${CMAKE_DEBUG_POSTFIX}")

target_link_libraries(${PROJECT_NAME} PUBLIC "mia")
target_link_libraries(${PROJECT_NAME} PUBLIC "itkio")
target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_PROGRAM_OPTIONS_LIBRARY})