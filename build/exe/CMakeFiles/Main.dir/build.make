# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/Frederik/Documents/C++/Projects/Interpolate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/Frederik/Documents/C++/Projects/Interpolate/build

# Include any dependencies generated for this target.
include exe/CMakeFiles/Main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include exe/CMakeFiles/Main.dir/compiler_depend.make

# Include the progress variables for this target.
include exe/CMakeFiles/Main.dir/progress.make

# Include the compile flags for this target's objects.
include exe/CMakeFiles/Main.dir/flags.make

exe/CMakeFiles/Main.dir/main.cpp.o: exe/CMakeFiles/Main.dir/flags.make
exe/CMakeFiles/Main.dir/main.cpp.o: /Users/Frederik/Documents/C++/Projects/Interpolate/exe/main.cpp
exe/CMakeFiles/Main.dir/main.cpp.o: exe/CMakeFiles/Main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/Frederik/Documents/C++/Projects/Interpolate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object exe/CMakeFiles/Main.dir/main.cpp.o"
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT exe/CMakeFiles/Main.dir/main.cpp.o -MF CMakeFiles/Main.dir/main.cpp.o.d -o CMakeFiles/Main.dir/main.cpp.o -c /Users/Frederik/Documents/C++/Projects/Interpolate/exe/main.cpp

exe/CMakeFiles/Main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Main.dir/main.cpp.i"
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/Frederik/Documents/C++/Projects/Interpolate/exe/main.cpp > CMakeFiles/Main.dir/main.cpp.i

exe/CMakeFiles/Main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Main.dir/main.cpp.s"
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/Frederik/Documents/C++/Projects/Interpolate/exe/main.cpp -o CMakeFiles/Main.dir/main.cpp.s

# Object files for target Main
Main_OBJECTS = \
"CMakeFiles/Main.dir/main.cpp.o"

# External object files for target Main
Main_EXTERNAL_OBJECTS =

/Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main: exe/CMakeFiles/Main.dir/main.cpp.o
/Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main: exe/CMakeFiles/Main.dir/build.make
/Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main: src/libinterpolateLib.a
/Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main: exe/CMakeFiles/Main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/Frederik/Documents/C++/Projects/Interpolate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main"
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
exe/CMakeFiles/Main.dir/build: /Users/Frederik/Documents/C++/Projects/Interpolate/exe/Main
.PHONY : exe/CMakeFiles/Main.dir/build

exe/CMakeFiles/Main.dir/clean:
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe && $(CMAKE_COMMAND) -P CMakeFiles/Main.dir/cmake_clean.cmake
.PHONY : exe/CMakeFiles/Main.dir/clean

exe/CMakeFiles/Main.dir/depend:
	cd /Users/Frederik/Documents/C++/Projects/Interpolate/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/Frederik/Documents/C++/Projects/Interpolate /Users/Frederik/Documents/C++/Projects/Interpolate/exe /Users/Frederik/Documents/C++/Projects/Interpolate/build /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe /Users/Frederik/Documents/C++/Projects/Interpolate/build/exe/CMakeFiles/Main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : exe/CMakeFiles/Main.dir/depend

