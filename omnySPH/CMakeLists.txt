cmake_minimum_required(VERSION 3.15...3.22)
cmake_policy(SET CMP0015 NEW)

project(omnySPH VERSION "0.0.2")

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

# Now we can find pybind11
find_package(pybind11 CONFIG REQUIRED)
# find_package(yaml-cpp REQUIRED)

pybind11_add_module(_core MODULE src/main.cpp)

link_directories("${CMAKE_SOURCE_DIR}/../build/lib")
message(WARNING "${CMAKE_SOURCE_DIR}/../build/lib")
target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
target_link_libraries(_core PRIVATE stdc++ m "${CMAKE_SOURCE_DIR}/../build/lib/libsimulation.a" "${CMAKE_SOURCE_DIR}/../build/lib/libtools.a")

install(TARGETS _core DESTINATION .)
