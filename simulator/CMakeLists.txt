CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(PathControlSimulator)

# Options
SET(EIGEN3_INCLUDE_DIR "Eigen 3 Include Directory." CACHE STRING "")

# Configure Compiler
ADD_COMPILE_OPTIONS(-std=c++11)

# Dependencies
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}" "include/")
FIND_PACKAGE(Boost 1.54 REQUIRED)

# Source Files
SET(SimSrc	src/main.cpp
			src/diffdrive.cpp
			src/spline.cpp
			src/quinticLevelSet.cpp
			src/quinticJacobian.cpp
			src/quinticHessian.cpp
			src/sfcontrol.cpp)

# Generate Executable
ADD_EXECUTABLE(followingsim ${SimSrc})
TARGET_INCLUDE_DIRECTORIES(followingsim INTERFACE include)
