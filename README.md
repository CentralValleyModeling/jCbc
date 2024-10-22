# Introduction

jCbc is a Java Native Interface for COIN OR Mixed Integer Linear Programming Solver Cbc and Linear Programming Solver CLP. jCbc  uses open source Simplified Wrapper and Interface Generator SWIG, which is a software development tool that connects libraries written in C and C++ with a variety of high-level programming languages. The MILP can be either constructed using the relavant functions in jCbc, or it can be loaded from an LP or MPS file. In either case, jCbc uses branch and cut algorithm to find an optimal solution or declares that the problem has no solution. 



*********************************************************************************
# License 

jCbc is distributed under GPL-v3 licese. Everyone is permitted to copy and distribute verbatim copies of this license document, but changing it is not allowed.

# Dependencies

The following dependencies are needed :-
* JDK version 1.8x
* Visual Studio 2022 + Intel Oneapi compilers
* SWIG version 4.2.1
*   Cbc release 2.9.x
*   Cbc release 2.10.x

Not required but recommended
* Cmake version 3.30+
* Git
* VSCode

# Compilation

 run_cmake.bat for details and CMakeLists.txt may need to be changed to your install locations and dependencies.

 Once run_cmake.bat is run it will create a solution file jCbc.sln in the build folder. Visual Studio can then use SWIG to generate the wrappers for C++ and Java and compile the code into jCbc.dll

  The generated files in build can then be compiled using a java compiler and packaged in a jar file (jCbc.jar)

  These two files (jCbc.dll and jCbc.jar) can be used to run the examples under examples directory and verify the results with values in *.lp files.
