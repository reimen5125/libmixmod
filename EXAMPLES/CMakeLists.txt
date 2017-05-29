# Configuration
#--------------
cmake_minimum_required(VERSION 2.8.5)
PROJECT(mixmod_EXAMPLES)


# Switch: cmake .. -DCMAKE_INSTALL_PREFIX="/usr/local/" ou "C:\Windows\ProgramFiles\mixmod\" ou ...
# NOTE: should be copied in every individual source package

# Installation prefix
#---------------------
if (NOT DEFINED CMAKE_INSTALL_PREFIX)
	set(CMAKE_INSTALL_PREFIX "/usr/local")
endif ()
message("CMAKE_INSTALL_PREFIX : " ${CMAKE_INSTALL_PREFIX})
# 
if (NOT DEFINED CPACK_PACKAGING_INSTALL_PREFIX)
	set(CPACK_PACKAGING_INSTALL_PREFIX "/usr")
endif ()
message("CPACK_PACKAGING_INSTALL_PREFIX : " ${CPACK_PACKAGING_INSTALL_PREFIX})

# Release/Debug
# #--------------
# if (NOT DEFINED CMAKE_BUILD_TYPE)
# 	set(CMAKE_BUILD_TYPE "Release")
# endif ()
# message(STATUS ${CMAKE_BUILD_TYPE})
IF(NOT  CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(CMAKE_BUILD_TYPE "Release")
endif ()
message("CMAKE_BUILD_TYPE : " ${CMAKE_BUILD_TYPE})

# Compiler settings
#-------------------
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g3")
#add the right flags for C++11 features [TODO: check if OK with MinGW/Windows]
# execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
#                 OUTPUT_VARIABLE GCC_VERSION)
# if (GCC_VERSION VERSION_GREATER 4.6)
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
# else ()
# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
# endif ()
message("CMAKE_CXX_FLAGS : " ${CMAKE_CXX_FLAGS})
#fix related to CMP0009. See cmake --help-policy CMP0009
#NOTE: we should not need it if required cmake VERSION >= 2.8.X... (?!)
if(POLICY CMP0009)
	cmake_policy(SET CMP0009 NEW)
endif(POLICY CMP0009)


# Headers and libraries
#--------------
#TODO
set(mixmodLib_include ${CMAKE_INSTALL_PREFIX}/include)
set(mixmodNEWMAT_include ${CMAKE_INSTALL_PREFIX}/include)

include_directories(${mixmodLib_include})
include_directories(${mixmodNEWMAT_include}/NEWMAT)

find_library(libmixmod mixmod PATHS ${CPACK_PACKAGING_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/lib NO_DEFAULT_PATH)
message("libmixmod found here : "  ${libmixmod})
find_library(libmixmod_newmat mixmod_newmat PATHS ${CPACK_PACKAGING_INSTALL_PREFIX}/lib ${CMAKE_INSTALL_PREFIX}/lib NO_DEFAULT_PATH)
message("libmixmod_newmat found here : "  ${libmixmod_newmat})

#sources
file(GLOB sources ${CMAKE_SOURCE_DIR}/main.cpp ${CMAKE_SOURCE_DIR}/clustering/*.cpp ${CMAKE_SOURCE_DIR}/discriminant_analysis/*.cpp)
add_executable(mixmod_example ${sources})
target_link_libraries(mixmod_example ${libmixmod} ${libmixmod_newmat})
#On install:
INSTALL(TARGETS mixmod_example DESTINATION bin COMPONENT dev)
