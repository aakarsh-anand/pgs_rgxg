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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev

# Include any dependencies generated for this target.
include CMakeFiles/GENIE_GxG.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/GENIE_GxG.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GENIE_GxG.dir/flags.make

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o: CMakeFiles/GENIE_GxG.dir/flags.make
CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o: ../src/gg_asymp_se_dev.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o"
	/u/local/compilers/gcc/7.5.0/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o -c /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/gg_asymp_se_dev.cpp

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.i"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/gg_asymp_se_dev.cpp > CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.i

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.s"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/gg_asymp_se_dev.cpp -o CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.s

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.requires:
.PHONY : CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.requires

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.provides: CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.requires
	$(MAKE) -f CMakeFiles/GENIE_GxG.dir/build.make CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.provides.build
.PHONY : CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.provides

CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.provides.build: CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o: CMakeFiles/GENIE_GxG.dir/flags.make
CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o: ../src/genotype.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o"
	/u/local/compilers/gcc/7.5.0/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o -c /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/genotype.cpp

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.i"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/genotype.cpp > CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.i

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.s"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/genotype.cpp -o CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.s

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.requires:
.PHONY : CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.requires

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.provides: CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.requires
	$(MAKE) -f CMakeFiles/GENIE_GxG.dir/build.make CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.provides.build
.PHONY : CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.provides

CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.provides.build: CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o: CMakeFiles/GENIE_GxG.dir/flags.make
CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o: ../src/storage.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o"
	/u/local/compilers/gcc/7.5.0/bin/g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o -c /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/storage.cpp

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GENIE_GxG.dir/src/storage.cpp.i"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/storage.cpp > CMakeFiles/GENIE_GxG.dir/src/storage.cpp.i

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GENIE_GxG.dir/src/storage.cpp.s"
	/u/local/compilers/gcc/7.5.0/bin/g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/src/storage.cpp -o CMakeFiles/GENIE_GxG.dir/src/storage.cpp.s

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.requires:
.PHONY : CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.requires

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.provides: CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.requires
	$(MAKE) -f CMakeFiles/GENIE_GxG.dir/build.make CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.provides.build
.PHONY : CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.provides

CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.provides.build: CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o

# Object files for target GENIE_GxG
GENIE_GxG_OBJECTS = \
"CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o" \
"CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o" \
"CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o"

# External object files for target GENIE_GxG
GENIE_GxG_EXTERNAL_OBJECTS =

GENIE_GxG: CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o
GENIE_GxG: CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o
GENIE_GxG: CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o
GENIE_GxG: CMakeFiles/GENIE_GxG.dir/build.make
GENIE_GxG: CMakeFiles/GENIE_GxG.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable GENIE_GxG"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GENIE_GxG.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GENIE_GxG.dir/build: GENIE_GxG
.PHONY : CMakeFiles/GENIE_GxG.dir/build

CMakeFiles/GENIE_GxG.dir/requires: CMakeFiles/GENIE_GxG.dir/src/gg_asymp_se_dev.cpp.o.requires
CMakeFiles/GENIE_GxG.dir/requires: CMakeFiles/GENIE_GxG.dir/src/genotype.cpp.o.requires
CMakeFiles/GENIE_GxG.dir/requires: CMakeFiles/GENIE_GxG.dir/src/storage.cpp.o.requires
.PHONY : CMakeFiles/GENIE_GxG.dir/requires

CMakeFiles/GENIE_GxG.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GENIE_GxG.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GENIE_GxG.dir/clean

CMakeFiles/GENIE_GxG.dir/depend:
	cd /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG_dev/GENIE_gxg/build_dev/CMakeFiles/GENIE_GxG.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GENIE_GxG.dir/depend

