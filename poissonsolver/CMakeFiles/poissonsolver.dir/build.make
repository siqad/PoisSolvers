# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/nathan/git/PoisSolvers/poissonsolver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nathan/git/PoisSolvers/poissonsolver

# Include any dependencies generated for this target.
include CMakeFiles/poissonsolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/poissonsolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/poissonsolver.dir/flags.make

CMakeFiles/poissonsolver.dir/src/solver.cpp.o: CMakeFiles/poissonsolver.dir/flags.make
CMakeFiles/poissonsolver.dir/src/solver.cpp.o: src/solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nathan/git/PoisSolvers/poissonsolver/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/poissonsolver.dir/src/solver.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonsolver.dir/src/solver.cpp.o -c /home/nathan/git/PoisSolvers/poissonsolver/src/solver.cpp

CMakeFiles/poissonsolver.dir/src/solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonsolver.dir/src/solver.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nathan/git/PoisSolvers/poissonsolver/src/solver.cpp > CMakeFiles/poissonsolver.dir/src/solver.cpp.i

CMakeFiles/poissonsolver.dir/src/solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonsolver.dir/src/solver.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nathan/git/PoisSolvers/poissonsolver/src/solver.cpp -o CMakeFiles/poissonsolver.dir/src/solver.cpp.s

CMakeFiles/poissonsolver.dir/src/solver.cpp.o.requires:

.PHONY : CMakeFiles/poissonsolver.dir/src/solver.cpp.o.requires

CMakeFiles/poissonsolver.dir/src/solver.cpp.o.provides: CMakeFiles/poissonsolver.dir/src/solver.cpp.o.requires
	$(MAKE) -f CMakeFiles/poissonsolver.dir/build.make CMakeFiles/poissonsolver.dir/src/solver.cpp.o.provides.build
.PHONY : CMakeFiles/poissonsolver.dir/src/solver.cpp.o.provides

CMakeFiles/poissonsolver.dir/src/solver.cpp.o.provides.build: CMakeFiles/poissonsolver.dir/src/solver.cpp.o


CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o: CMakeFiles/poissonsolver.dir/flags.make
CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o: src/electrodes.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nathan/git/PoisSolvers/poissonsolver/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o -c /home/nathan/git/PoisSolvers/poissonsolver/src/electrodes.cpp

CMakeFiles/poissonsolver.dir/src/electrodes.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonsolver.dir/src/electrodes.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nathan/git/PoisSolvers/poissonsolver/src/electrodes.cpp > CMakeFiles/poissonsolver.dir/src/electrodes.cpp.i

CMakeFiles/poissonsolver.dir/src/electrodes.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonsolver.dir/src/electrodes.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nathan/git/PoisSolvers/poissonsolver/src/electrodes.cpp -o CMakeFiles/poissonsolver.dir/src/electrodes.cpp.s

CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.requires:

.PHONY : CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.requires

CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.provides: CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.requires
	$(MAKE) -f CMakeFiles/poissonsolver.dir/build.make CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.provides.build
.PHONY : CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.provides

CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.provides.build: CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o


CMakeFiles/poissonsolver.dir/main.cpp.o: CMakeFiles/poissonsolver.dir/flags.make
CMakeFiles/poissonsolver.dir/main.cpp.o: main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nathan/git/PoisSolvers/poissonsolver/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/poissonsolver.dir/main.cpp.o"
	/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/poissonsolver.dir/main.cpp.o -c /home/nathan/git/PoisSolvers/poissonsolver/main.cpp

CMakeFiles/poissonsolver.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/poissonsolver.dir/main.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nathan/git/PoisSolvers/poissonsolver/main.cpp > CMakeFiles/poissonsolver.dir/main.cpp.i

CMakeFiles/poissonsolver.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/poissonsolver.dir/main.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nathan/git/PoisSolvers/poissonsolver/main.cpp -o CMakeFiles/poissonsolver.dir/main.cpp.s

CMakeFiles/poissonsolver.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/poissonsolver.dir/main.cpp.o.requires

CMakeFiles/poissonsolver.dir/main.cpp.o.provides: CMakeFiles/poissonsolver.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/poissonsolver.dir/build.make CMakeFiles/poissonsolver.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/poissonsolver.dir/main.cpp.o.provides

CMakeFiles/poissonsolver.dir/main.cpp.o.provides.build: CMakeFiles/poissonsolver.dir/main.cpp.o


# Object files for target poissonsolver
poissonsolver_OBJECTS = \
"CMakeFiles/poissonsolver.dir/src/solver.cpp.o" \
"CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o" \
"CMakeFiles/poissonsolver.dir/main.cpp.o"

# External object files for target poissonsolver
poissonsolver_EXTERNAL_OBJECTS =

poissonsolver: CMakeFiles/poissonsolver.dir/src/solver.cpp.o
poissonsolver: CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o
poissonsolver: CMakeFiles/poissonsolver.dir/main.cpp.o
poissonsolver: CMakeFiles/poissonsolver.dir/build.make
poissonsolver: CMakeFiles/poissonsolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nathan/git/PoisSolvers/poissonsolver/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable poissonsolver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/poissonsolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/poissonsolver.dir/build: poissonsolver

.PHONY : CMakeFiles/poissonsolver.dir/build

CMakeFiles/poissonsolver.dir/requires: CMakeFiles/poissonsolver.dir/src/solver.cpp.o.requires
CMakeFiles/poissonsolver.dir/requires: CMakeFiles/poissonsolver.dir/src/electrodes.cpp.o.requires
CMakeFiles/poissonsolver.dir/requires: CMakeFiles/poissonsolver.dir/main.cpp.o.requires

.PHONY : CMakeFiles/poissonsolver.dir/requires

CMakeFiles/poissonsolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/poissonsolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/poissonsolver.dir/clean

CMakeFiles/poissonsolver.dir/depend:
	cd /home/nathan/git/PoisSolvers/poissonsolver && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nathan/git/PoisSolvers/poissonsolver /home/nathan/git/PoisSolvers/poissonsolver /home/nathan/git/PoisSolvers/poissonsolver /home/nathan/git/PoisSolvers/poissonsolver /home/nathan/git/PoisSolvers/poissonsolver/CMakeFiles/poissonsolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/poissonsolver.dir/depend

