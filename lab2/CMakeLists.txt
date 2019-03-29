cmake_minimum_required (VERSION 2.6)
project (Mownit_lab2)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_library (Matrix src/matrix.c)
add_library (Matrixcpp src/matrixcpp.cpp)
add_library (Circuit src/circuit.cpp)


# add the executable
add_executable (zad1 src/zad1.c)
add_executable (zad2 src/zad2.c)
add_executable (  zad3 src/zad3.cpp)


target_link_libraries (zad1 Matrix)
target_link_libraries (zad2 Matrix)
target_link_libraries (zad3 Matrixcpp Circuit)
