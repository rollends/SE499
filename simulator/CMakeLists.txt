
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(PathControlSimulator)

# Options
SET(EIGEN3_INCLUDE_DIR "Eigen 3 Include Directory." CACHE STRING "")
SET(BOOST_1_63_0_INCLUDE_DIR "Boost 1.63.0 Header-Only Include Directory." CACHE STRING "")

# Configure Compiler
ADD_COMPILE_OPTIONS(-std=c++11)

# Dependencies
INCLUDE_DIRECTORIES( "${EIGEN3_INCLUDE_DIR}"
                     "${BOOST_1_63_0_INCLUDE_DIR}"
                     "include/" )

# Source Files
SET(SimSrc  src/diffdrive.cpp
            src/spline.cpp
            src/quinticLevelSet.cpp
            src/quinticJacobian.cpp
            src/quinticHessian.cpp
            src/sfcontrol.cpp
            src/world.cpp
            src/plan.cpp)

SET(TestFiles   src/tests/testSpline.cpp
                src/tests/testSerretFrenetControl.cpp
                src/tests/testWorld.cpp)

# Generate Temporary Static lib
ADD_LIBRARY(rose499 STATIC ${SimSrc})

# Generate Executable
ADD_EXECUTABLE(followingsim src/main.cpp)
TARGET_LINK_LIBRARIES(followingsim rose499)
ADD_EXECUTABLE(testfollowingsim src/tests/testMain.cpp ${TestFiles})
TARGET_LINK_LIBRARIES(testfollowingsim rose499)
