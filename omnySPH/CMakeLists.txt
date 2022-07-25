cmake_minimum_required(VERSION 3.15...3.22)
cmake_policy(SET CMP0015 NEW)

project(omnySPH VERSION "0.0.2")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/cmake")
# SET(CMAKE_PREFIX_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR};${CMAKE_SOURCE_DIR}/cmake")

SET(EIGEN3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../vcpkg/installed/x64-linux/include")
# SET(YAML_CPP_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/../../vcpkg/installed/x64-linux/include")

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific directory
  # here.
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp -DNOMINMAX -D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -fpic")

SET(INCLUDE_DIRS
	${CMAKE_SOURCE_DIR} 	
	"${CMAKE_SOURCE_DIR}/../"
	${CMAKE_SOURCE_DIR} 		
	${CMAKE_BINARY_DIR}   
	# ${YAML_CPP_INCLUDE_DIR}
	# ${GLAD_INCLUDE_DIR}
	# ${GLFW3_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
	# ${BOOST_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})   

# Now we can find pybind11
find_package(pybind11 CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)
# find_package(yaml-cpp REQUIRED)

pybind11_add_module(_core MODULE src/main.cpp)

link_directories("${CMAKE_SOURCE_DIR}/../build/lib")
message(WARNING "${CMAKE_SOURCE_DIR}/../build/lib")
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_core PRIVATE stdc++ m 
  "${CMAKE_SOURCE_DIR}/../build/lib/libsimulation.a" "${CMAKE_SOURCE_DIR}/../build/lib/libtools.a" "${CMAKE_SOURCE_DIR}/../build/lib/libimgui.a" 
  "${CMAKE_SOURCE_DIR}/../../vcpkg/installed/x64-linux/lib/libyaml-cpp.a"
  )

install(TARGETS _core DESTINATION .)
