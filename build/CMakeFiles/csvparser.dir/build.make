# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build

# Include any dependencies generated for this target.
include CMakeFiles/csvparser.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/csvparser.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/csvparser.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/csvparser.dir/flags.make

CMakeFiles/csvparser.dir/src/csvparser.cpp.o: CMakeFiles/csvparser.dir/flags.make
CMakeFiles/csvparser.dir/src/csvparser.cpp.o: /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csvparser.cpp
CMakeFiles/csvparser.dir/src/csvparser.cpp.o: CMakeFiles/csvparser.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/csvparser.dir/src/csvparser.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/csvparser.dir/src/csvparser.cpp.o -MF CMakeFiles/csvparser.dir/src/csvparser.cpp.o.d -o CMakeFiles/csvparser.dir/src/csvparser.cpp.o -c /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csvparser.cpp

CMakeFiles/csvparser.dir/src/csvparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/csvparser.dir/src/csvparser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csvparser.cpp > CMakeFiles/csvparser.dir/src/csvparser.cpp.i

CMakeFiles/csvparser.dir/src/csvparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/csvparser.dir/src/csvparser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csvparser.cpp -o CMakeFiles/csvparser.dir/src/csvparser.cpp.s

CMakeFiles/csvparser.dir/src/csv_bind.cpp.o: CMakeFiles/csvparser.dir/flags.make
CMakeFiles/csvparser.dir/src/csv_bind.cpp.o: /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csv_bind.cpp
CMakeFiles/csvparser.dir/src/csv_bind.cpp.o: CMakeFiles/csvparser.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/csvparser.dir/src/csv_bind.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/csvparser.dir/src/csv_bind.cpp.o -MF CMakeFiles/csvparser.dir/src/csv_bind.cpp.o.d -o CMakeFiles/csvparser.dir/src/csv_bind.cpp.o -c /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csv_bind.cpp

CMakeFiles/csvparser.dir/src/csv_bind.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/csvparser.dir/src/csv_bind.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csv_bind.cpp > CMakeFiles/csvparser.dir/src/csv_bind.cpp.i

CMakeFiles/csvparser.dir/src/csv_bind.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/csvparser.dir/src/csv_bind.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/src/csv_bind.cpp -o CMakeFiles/csvparser.dir/src/csv_bind.cpp.s

# Object files for target csvparser
csvparser_OBJECTS = \
"CMakeFiles/csvparser.dir/src/csvparser.cpp.o" \
"CMakeFiles/csvparser.dir/src/csv_bind.cpp.o"

# External object files for target csvparser
csvparser_EXTERNAL_OBJECTS =

csvparser.cpython-311-x86_64-linux-gnu.so: CMakeFiles/csvparser.dir/src/csvparser.cpp.o
csvparser.cpython-311-x86_64-linux-gnu.so: CMakeFiles/csvparser.dir/src/csv_bind.cpp.o
csvparser.cpython-311-x86_64-linux-gnu.so: CMakeFiles/csvparser.dir/build.make
csvparser.cpython-311-x86_64-linux-gnu.so: CMakeFiles/csvparser.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared module csvparser.cpython-311-x86_64-linux-gnu.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/csvparser.dir/link.txt --verbose=$(VERBOSE)
	/usr/bin/strip /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build/csvparser.cpython-311-x86_64-linux-gnu.so

# Rule to build all files generated by this target.
CMakeFiles/csvparser.dir/build: csvparser.cpython-311-x86_64-linux-gnu.so
.PHONY : CMakeFiles/csvparser.dir/build

CMakeFiles/csvparser.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/csvparser.dir/cmake_clean.cmake
.PHONY : CMakeFiles/csvparser.dir/clean

CMakeFiles/csvparser.dir/depend:
	cd /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3 /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3 /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build /home/andreg/Desktop/uni/ap/homeworks/hm3/AP_homework3/build/CMakeFiles/csvparser.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/csvparser.dir/depend
