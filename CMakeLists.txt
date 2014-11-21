cmake_minimum_required(VERSION 2.8.3)
project(aslam_sample_consensus)

find_package(catkin_simple REQUIRED)

catkin_simple()

# Helen's block of happy fun make this stuff run.
if(NOT DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_definitions(-O3)
add_definitions(-w)

# Avoid clash with tr1::tuple: https://code.google.com/p/googletest/source/browse/trunk/README?r=589#257
add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)

#catkin_add_gtest(${PROJECT_NAME}_test
#  test/test_main.cpp
#  test/TestRansac.cpp
#)

#TARGET_LINK_LIBRARIES( ${PROJECT_NAME}_test
#pthread
#)
cs_install()
cs_export()
