# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /var/lib/snapd/snap/clion/162/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /var/lib/snapd/snap/clion/162/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tom/Escritorio/discrete-optimization/tsp/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/tsp.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/tsp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tsp.dir/flags.make

CMakeFiles/tsp.dir/main.cpp.o: CMakeFiles/tsp.dir/flags.make
CMakeFiles/tsp.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tsp.dir/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tsp.dir/main.cpp.o -c /home/tom/Escritorio/discrete-optimization/tsp/src/main.cpp

CMakeFiles/tsp.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tsp.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tom/Escritorio/discrete-optimization/tsp/src/main.cpp > CMakeFiles/tsp.dir/main.cpp.i

CMakeFiles/tsp.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tsp.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tom/Escritorio/discrete-optimization/tsp/src/main.cpp -o CMakeFiles/tsp.dir/main.cpp.s

CMakeFiles/tsp.dir/Graph/Graph.cpp.o: CMakeFiles/tsp.dir/flags.make
CMakeFiles/tsp.dir/Graph/Graph.cpp.o: ../Graph/Graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/tsp.dir/Graph/Graph.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tsp.dir/Graph/Graph.cpp.o -c /home/tom/Escritorio/discrete-optimization/tsp/src/Graph/Graph.cpp

CMakeFiles/tsp.dir/Graph/Graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tsp.dir/Graph/Graph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tom/Escritorio/discrete-optimization/tsp/src/Graph/Graph.cpp > CMakeFiles/tsp.dir/Graph/Graph.cpp.i

CMakeFiles/tsp.dir/Graph/Graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tsp.dir/Graph/Graph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tom/Escritorio/discrete-optimization/tsp/src/Graph/Graph.cpp -o CMakeFiles/tsp.dir/Graph/Graph.cpp.s

CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o: CMakeFiles/tsp.dir/flags.make
CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o: ../UnionFind/UnionFind.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o -c /home/tom/Escritorio/discrete-optimization/tsp/src/UnionFind/UnionFind.cpp

CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tom/Escritorio/discrete-optimization/tsp/src/UnionFind/UnionFind.cpp > CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.i

CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tom/Escritorio/discrete-optimization/tsp/src/UnionFind/UnionFind.cpp -o CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.s

CMakeFiles/tsp.dir/Utils/utils.cpp.o: CMakeFiles/tsp.dir/flags.make
CMakeFiles/tsp.dir/Utils/utils.cpp.o: ../Utils/utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/tsp.dir/Utils/utils.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tsp.dir/Utils/utils.cpp.o -c /home/tom/Escritorio/discrete-optimization/tsp/src/Utils/utils.cpp

CMakeFiles/tsp.dir/Utils/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tsp.dir/Utils/utils.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/tom/Escritorio/discrete-optimization/tsp/src/Utils/utils.cpp > CMakeFiles/tsp.dir/Utils/utils.cpp.i

CMakeFiles/tsp.dir/Utils/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tsp.dir/Utils/utils.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/tom/Escritorio/discrete-optimization/tsp/src/Utils/utils.cpp -o CMakeFiles/tsp.dir/Utils/utils.cpp.s

# Object files for target tsp
tsp_OBJECTS = \
"CMakeFiles/tsp.dir/main.cpp.o" \
"CMakeFiles/tsp.dir/Graph/Graph.cpp.o" \
"CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o" \
"CMakeFiles/tsp.dir/Utils/utils.cpp.o"

# External object files for target tsp
tsp_EXTERNAL_OBJECTS =

tsp: CMakeFiles/tsp.dir/main.cpp.o
tsp: CMakeFiles/tsp.dir/Graph/Graph.cpp.o
tsp: CMakeFiles/tsp.dir/UnionFind/UnionFind.cpp.o
tsp: CMakeFiles/tsp.dir/Utils/utils.cpp.o
tsp: CMakeFiles/tsp.dir/build.make
tsp: CMakeFiles/tsp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable tsp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tsp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tsp.dir/build: tsp
.PHONY : CMakeFiles/tsp.dir/build

CMakeFiles/tsp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tsp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tsp.dir/clean

CMakeFiles/tsp.dir/depend:
	cd /home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tom/Escritorio/discrete-optimization/tsp/src /home/tom/Escritorio/discrete-optimization/tsp/src /home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug /home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug /home/tom/Escritorio/discrete-optimization/tsp/src/cmake-build-debug/CMakeFiles/tsp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tsp.dir/depend

