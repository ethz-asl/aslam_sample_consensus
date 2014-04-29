cmake_minimum_required(VERSION 2.8.3)
project(aslam_sample_consensus)

find_package(catkin_simple REQUIRED)

catkin_simple()

find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIRS})
add_definitions(${EIGEN_DEFINITIONS})

catkin_add_gtest(${PROJECT_NAME}_test
  test/test_main.cpp
  test/TestRansac.cpp
)

cs_install()
cs_export()
