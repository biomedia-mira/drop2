set(THIRD_PARTY_DIR $ENV{THIRD_PARTY_DIR})

list(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
list(APPEND CMAKE_PREFIX_PATH "${THIRD_PARTY_DIR}/qt5")
list(REMOVE_DUPLICATES CMAKE_PREFIX_PATH )

# the ITK library directory which contains the ITKConfig.cmake file
file(GLOB_RECURSE _USE_ITK_CMAKE ${THIRD_PARTY_DIR}/itk/lib/cmake ${THIRD_PARTY_DIR}/itk/lib/cmake/*/ITKConfig.cmake)
list(LENGTH _USE_ITK_CMAKE _NUM_USE_ITK_CMAKE_FILES)
if (${_NUM_USE_ITK_CMAKE_FILES} EQUAL 1)
  get_filename_component(ITK_DIR ${_USE_ITK_CMAKE} PATH) 
endif()

if(WIN32)

	# Boost
	set(Boost_USE_STATIC_LIBS ON)
	set(Boost_USE_MULTITHREADED ON)
	set(Boost_USE_STATIC_RUNTIME OFF)
	set(BOOST_ROOT ${THIRD_PARTY_DIR}/boost)

	# OpenCV
	set(OpenCV_DIR "${THIRD_PARTY_DIR}/opencv")

	# TBB
	set(TBB_TARGET_VS "")
	if(MSVC11)
    set(TBB_TARGET_VS "vc11")
	elseif(MSVC12)
    set(TBB_TARGET_VS "vc12")
	elseif(MSVC14)
    set(TBB_TARGET_VS "vc14")
	endif()

	set(ENV{TBB_ARCH_PLATFORM} "intel64/${TBB_TARGET_VS}")
	set(TBB_INSTALL_DIR "${THIRD_PARTY_DIR}/tbb")

	set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${THIRD_PARTY_DIR}/glew/lib)

	# zlib/minizip
	set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${THIRD_PARTY_DIR}/zlib/include)
	set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${THIRD_PARTY_DIR}/zlib/lib)

	set(TIFF_NAMES tiff libtiff itktiff)
	set(JPEG_NAMES jpeg libjpeg itkjpeg itkjpeg8 itkjpeg12 itkjpeg16)
	set(ZLIB_NAMES zlib itkzlib)
	set(OPENJPEG_NAMES openjpeg itkopenjpeg)

	# hdf5
	list(APPEND CMAKE_MODULE_PATH "${THIRD_PARTY_DIR}/hdf5/cmake/hdf5")
	list(APPEND CMAKE_FRAMEWORK_PATH ${THIRD_PARTY_DIR}/hdf5)
endif()