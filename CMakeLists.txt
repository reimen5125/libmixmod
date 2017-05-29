# Configuration
#--------------
cmake_minimum_required(VERSION 2.8.5)
PROJECT(mixmod)

#usgae "cmake . " + Options :
#	-DCMAKE_INSTALL_PREFIX="<Path>" 
#		Default : /usr/local
#		-> /usr/local/lib : libmixmod.so, libmixmod_newmat.so (module lib)
#		-> /usr/local/include/mixmod : *.h  + /usr/local/share/mixmod : EXAMPLES (module dev)
#	-DCMAKE_BUILD_TYPE="Debug" Default="Release"

# Switch: cmake .. -DCMAKE_INSTALL_PREFIX="/usr/" ou "C:\Windows\ProgramFiles\mixmod\" ou ...
# NOTE: should be copied in every individual source package

# cpack -G DEB

# version mumber
set(CPACK_PACKAGE_VERSION_MAJOR "3")
set(CPACK_PACKAGE_VERSION_MINOR "2")
set(CPACK_PACKAGE_VERSION_PATCH "2")

# Installation prefix, packaging prefix
#--------------------------------------
if (NOT DEFINED CMAKE_INSTALL_PREFIX)
	set(CMAKE_INSTALL_PREFIX "/usr/local")
endif ()
message("CMAKE_INSTALL_PREFIX :"  ${CMAKE_INSTALL_PREFIX})

if (NOT DEFINED CPACK_PACKAGING_INSTALL_PREFIX)
	set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
endif ()
message("CPACK_PACKAGING_INSTALL_PREFIX : " ${CPACK_PACKAGING_INSTALL_PREFIX})

# Release/Debug
#--------------
IF(NOT  CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(CMAKE_BUILD_TYPE "Release")
endif ()
message("CMAKE_BUILD_TYPE : " ${CMAKE_BUILD_TYPE})

# Compiler settings
#-------------------
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -Wno-deprecated")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wno-deprecated")
#find_package(OpenMP)
#if (OPENMP_FOUND)
#  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
#  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OpenMP_CXX_FLAGS}")
#endif()
message("CMAKE_CXX_FLAGS_RELEASE : " ${CMAKE_CXX_FLAGS_RELEASE})
message("CMAKE_CXX_FLAGS_DEBUG : " ${CMAKE_CXX_FLAGS_DEBUG})
message("CMAKE_CXX_FLAGS : " ${CMAKE_CXX_FLAGS})

#fix related to CMP0009. See cmake --help-policy CMP0009
#NOTE: we should not need it if required cmake VERSION >= 2.8.X... (?!)
if(POLICY CMP0009)
	cmake_policy(SET CMP0009 NEW)
endif(POLICY CMP0009)


#------------------------
# Recurse in directories
#------------------------
ADD_SUBDIRECTORY(SRC)
#ADD_SUBDIRECTORY(EXAMPLES)

#install some extra files
install(FILES "README" "LICENSE" 
	DESTINATION share/mixmod
	COMPONENT dev)

#-----
#CPack
#-----
include(${CMAKE_SOURCE_DIR}/CPackLists.txt)
include(CPack)
