cmake_minimum_required(VERSION 2.8.3)
project(aslam_sample_consensus)

find_package(catkin_simple REQUIRED)

catkin_simple()

find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
add_definitions(${EIGEN_DEFINITIONS})


# Avoid clash with tr1::tuple: https://code.google.com/p/googletest/source/browse/trunk/README?r=589#257
add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)

catkin_add_gtest(${PROJECT_NAME}_test
  test/test_main.cpp
  test/TestRansac.cpp
)

TARGET_LINK_LIBRARIES( ${PROJECT_NAME}_test
pthread
)
cs_install()
cs_export()
