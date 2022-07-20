cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
# CMake policy settings
cmake_policy(VERSION 3.9)												# required for boost 1.65.1
cmake_policy(SET CMP0071 NEW)											# automoc setting policy introduced in CMake 3.10
   
# CMake settings to build cuda projects properly and to work well with other depdencies, e.g. boost
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/cmake") # add local CMake files 
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON) 						
SET(CMAKE_INCLUDE_CURRENT_DIR ON)	# needed for QT 
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${CMAKE_BINARY_DIR}/bin)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/lib) 
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/bin)
     
# Project settings
SET(NAME omniSPH)
project(${NAME} LANGUAGES CXX) 

FIND_PACKAGE( yaml-cpp		REQUIRED )
FIND_PACKAGE( Eigen3			REQUIRED ) 

configure_file(config.h.in config/config.h) 
# Set Include paths, due to everything needing almost everything this is done here
SET(INCLUDE_DIRS
	${CMAKE_SOURCE_DIR} 	
	${CMAKE_SOURCE_DIR} 		
	${CMAKE_BINARY_DIR}   
	${YAML_CPP_INCLUDE_DIR}
	${GLAD_INCLUDE_DIR}
	${GLFW3_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})   
  
SET(GUI TRUE CACHE BOOL "User defined")
IF(NOT(GUI))
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp -DNO_GUI")
ELSE()
FIND_PACKAGE( glfw3			REQUIRED )
FIND_PACKAGE( glad			REQUIRED ) 
ENDIF()


IF(WIN32)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++latest /openmp /fpermissive- /arch:AVX2")
add_compile_definitions(NOMINMAX _USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
ELSE()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -DNOMINMAX -D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS")
ENDIF()



macro(_GLOB_SOURCES LIBRARY_NAME)
	file(GLOB_RECURSE CPP_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.cpp)
	file(GLOB_RECURSE CPP_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${LIBRARY_NAME}/*.h ${LIBRARY_NAME}/*.hpp)
	SET(source_list ${CPP_SOURCES} ${CPP_HEADERS})
	   
	foreach(source IN LISTS source_list)
		get_filename_component(source_path "${source}" PATH) 
		string(REPLACE "/" "\\" source_path_msvc "${source_path}")
		source_group("${source_path_msvc}" FILES "${source}")
	endforeach()  
	SET(source_list ${CPP_SOURCES} ${CPP_HEADERS})
endmacro(_GLOB_SOURCES) 
     
# Macro used to link all system libraries, e.g.  boost and TBB but NOT cuda
macro(_SETUP_DEPENDENCIES MODULE_NAME)
	target_link_libraries(${MODULE_NAME}
		yaml-cpp
		${ARGN}
	)  
	IF(WIN32)
	ELSE()	
		target_link_libraries(${MODULE_NAME} stdc++ m)
	ENDIF()
	target_include_directories(${MODULE_NAME} PUBLIC ${INCLUDE_DIRECTORIES})
endmacro(_SETUP_DEPENDENCIES)  
  
function(buildLibrary LIBRARY_NAME TARGETLIST PREFIX ) 
	_GLOB_SOURCES("${PREFIX}/${LIBRARY_NAME}")
	add_library( ${LIBRARY_NAME} ${CPP_SOURCES} ${CPP_HEADERS}) 
	_SETUP_DEPENDENCIES(${LIBRARY_NAME} ${ARGN})
endfunction(buildLibrary) 

IF(GUI)
_GLOB_SOURCES("gui")
SET(guiSources "${source_list}")
_GLOB_SOURCES("imgui")
SET(imguiSources "${source_list}")
ENDIF()
_GLOB_SOURCES("simulation")
SET(simSources "${source_list}")
_GLOB_SOURCES("tools")
SET(toolSources "${source_list}")

add_library(simulation ${simSources})
add_library(tools ${toolSources})
IF(GUI)
add_library(imgui ${imguiSources})

add_executable(gui Source.cpp ${guiSources})
_SETUP_DEPENDENCIES(imgui glfw glad::glad)
_SETUP_DEPENDENCIES(gui imgui glfw glad::glad simulation tools Eigen3::Eigen) 

add_executable(console SourceTerminal.cpp)

_SETUP_DEPENDENCIES(tools Eigen3::Eigen)
_SETUP_DEPENDENCIES(simulation tools Eigen3::Eigen)

_SETUP_DEPENDENCIES(console simulation tools Eigen3::Eigen imgui glfw glad::glad)
ELSE()
add_executable(console SourceTerminal.cpp)

_SETUP_DEPENDENCIES(tools Eigen3::Eigen)
_SETUP_DEPENDENCIES(simulation tools Eigen3::Eigen)

_SETUP_DEPENDENCIES(console simulation tools Eigen3::Eigen)

ENDIF()



# target_compile_options(rayTracingCore PRIVATE /arch:AVX2) 
