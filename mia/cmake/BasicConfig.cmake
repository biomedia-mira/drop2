cmake_policy(SET CMP0022 NEW)

#------------------------------------------------------------------------------
# Remove all configurations besides Debug and Release and set the debug
# postfix to "d".
#
if(WIN32)
	set(CMAKE_DEBUG_POSTFIX "d")
endif()

#------------------------------------------------------------------------------
# Build all binary files (.exe and .dll) into a single /bin folder.
#
if(WIN32)
  # Prevent accidental override of the output directory for submodules.
  if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
  endif()
endif()

#------------------------------------------------------------------------------
# Enable target grouping for MSVC
#
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

#------------------------------------------------------------------------------
# Configure an output postfix depending on the build type.
# TODO: verify that this is still required.
#
if(CMAKE_CL_64 OR CMAKE_GENERATOR MATCHES Win64)
else()
  add_definitions(-DEIGEN_DONT_ALIGN_STATICALLY)
endif()

#------------------------------------------------------------------------------
# This file is required to improve CMake's location capabilities.
#
include(FindConfig)

#------------------------------------------------------------------------------
# This file is required to improve the integration of 3rdparty libraries.
#
if (MSVC)
  # disable persistent and stupid warnings
  #
  # 4503: decorated name length exceeded, name was truncated
  # 4127: conditional expression is constant
  # 4251: class needs to have dll-interface to be used by clients of class
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4503 /wd4127 /wd4251")
  
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_CRT_SECURE_NO_WARNINGS /D_SCL_SECURE_NO_WARNINGS")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  
  if(CMAKE_CL_64)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} /bigobj")
  endif()
  
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS)
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_DEBUG)
  list(REMOVE_DUPLICATES CMAKE_CXX_FLAGS_RELWITHDEBINFO)

  add_definitions("-DNOMINMAX")
  add_definitions("-D_VARIADIC_MAX=10") # VC11 bug regarding std::pair
endif()

#------------------------------------------------------------------------------
# enable handling of structured exceptions
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHa")
endif(MSVC)

#------------------------------------------------------------------------------
# enabling C++ 11 support for gcc
if(CMAKE_COMPILER_IS_GNUCXX)
  list( APPEND CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} -O3")
endif()
