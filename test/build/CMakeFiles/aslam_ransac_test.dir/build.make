# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/2.8.8/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/2.8.8/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/local/Cellar/cmake/2.8.8/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build

# Include any dependencies generated for this target.
include CMakeFiles/aslam_ransac_test.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/aslam_ransac_test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/aslam_ransac_test.dir/flags.make

CMakeFiles/aslam_ransac_test.dir/test/test_main.o: CMakeFiles/aslam_ransac_test.dir/flags.make
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: ../test/test_main.cpp
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: ../manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_console/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_common/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_eigen/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: /Users/ptf/Work/ethz/code/ros_vc/eigen/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/test_main.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_random/manifest.xml
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/aslam_ransac_test.dir/test/test_main.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/aslam_ransac_test.dir/test/test_main.o -c /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/test_main.cpp

CMakeFiles/aslam_ransac_test.dir/test/test_main.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/aslam_ransac_test.dir/test/test_main.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/test_main.cpp > CMakeFiles/aslam_ransac_test.dir/test/test_main.i

CMakeFiles/aslam_ransac_test.dir/test/test_main.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/aslam_ransac_test.dir/test/test_main.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/test_main.cpp -o CMakeFiles/aslam_ransac_test.dir/test/test_main.s

CMakeFiles/aslam_ransac_test.dir/test/test_main.o.requires:
.PHONY : CMakeFiles/aslam_ransac_test.dir/test/test_main.o.requires

CMakeFiles/aslam_ransac_test.dir/test/test_main.o.provides: CMakeFiles/aslam_ransac_test.dir/test/test_main.o.requires
	$(MAKE) -f CMakeFiles/aslam_ransac_test.dir/build.make CMakeFiles/aslam_ransac_test.dir/test/test_main.o.provides.build
.PHONY : CMakeFiles/aslam_ransac_test.dir/test/test_main.o.provides

CMakeFiles/aslam_ransac_test.dir/test/test_main.o.provides.build: CMakeFiles/aslam_ransac_test.dir/test/test_main.o

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: CMakeFiles/aslam_ransac_test.dir/flags.make
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: ../test/TestRansac.cpp
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: ../manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_console/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_common/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_eigen/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: /Users/ptf/Work/ethz/code/ros_vc/eigen/manifest.xml
CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o: /Users/ptf/Work/ethz/code/ros_vc/schweizer_messer/sm_random/manifest.xml
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -o CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o -c /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/TestRansac.cpp

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/aslam_ransac_test.dir/test/TestRansac.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -E /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/TestRansac.cpp > CMakeFiles/aslam_ransac_test.dir/test/TestRansac.i

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/aslam_ransac_test.dir/test/TestRansac.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -W -Wall -Wno-unused-parameter -fno-strict-aliasing -pthread -S /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/test/TestRansac.cpp -o CMakeFiles/aslam_ransac_test.dir/test/TestRansac.s

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.requires:
.PHONY : CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.requires

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.provides: CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.requires
	$(MAKE) -f CMakeFiles/aslam_ransac_test.dir/build.make CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.provides.build
.PHONY : CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.provides

CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.provides.build: CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o

# Object files for target aslam_ransac_test
aslam_ransac_test_OBJECTS = \
"CMakeFiles/aslam_ransac_test.dir/test/test_main.o" \
"CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o"

# External object files for target aslam_ransac_test
aslam_ransac_test_EXTERNAL_OBJECTS =

../bin/aslam_ransac_test: CMakeFiles/aslam_ransac_test.dir/test/test_main.o
../bin/aslam_ransac_test: CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o
../bin/aslam_ransac_test: CMakeFiles/aslam_ransac_test.dir/build.make
../bin/aslam_ransac_test: CMakeFiles/aslam_ransac_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../bin/aslam_ransac_test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/aslam_ransac_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/aslam_ransac_test.dir/build: ../bin/aslam_ransac_test
.PHONY : CMakeFiles/aslam_ransac_test.dir/build

CMakeFiles/aslam_ransac_test.dir/requires: CMakeFiles/aslam_ransac_test.dir/test/test_main.o.requires
CMakeFiles/aslam_ransac_test.dir/requires: CMakeFiles/aslam_ransac_test.dir/test/TestRansac.o.requires
.PHONY : CMakeFiles/aslam_ransac_test.dir/requires

CMakeFiles/aslam_ransac_test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/aslam_ransac_test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/aslam_ransac_test.dir/clean

CMakeFiles/aslam_ransac_test.dir/depend:
	cd /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build /Users/ptf/Work/ethz/code/ros_vc/aslam_cv/aslam_ransac/build/CMakeFiles/aslam_ransac_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/aslam_ransac_test.dir/depend
