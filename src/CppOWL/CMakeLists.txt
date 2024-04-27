cmake_minimum_required(VERSION 3.10)
project(OWL_Ball_Toolbox)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Specify Boost root and find Boost
set(BOOST_ROOT "/usr/local/opt/boost")
find_package(Boost 1.84.0 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Include JlCxx package (assuming it is already installed); next line is overwritten with `CxxWrap.prefix()`
list(APPEND CMAKE_PREFIX_PATH "/Users/jakeroth/.julia/artifacts/2ce0502d68a1c277361b766c5050582bf167ad43")
find_package(JlCxx REQUIRED)

# Add your source files
set(SOURCE_FILES
    project_to_OWL_ball.cpp
    julia_interface.cpp
)

# Create a shared library
add_library(owlball SHARED ${SOURCE_FILES})

# Link JlCxx library
target_link_libraries(owlball JlCxx::cxxwrap_julia ${Boost_LIBRARIES})

# Specify include directories
target_include_directories(owlball PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})