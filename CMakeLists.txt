cmake_minimum_required(VERSION 3.3)
project(filterSqp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES test_adept.cpp)
add_executable(filterSqp ${SOURCE_FILES})
target_link_libraries(filterSqp adept)

add_executable(test_gsl test_gsl.cpp)
target_link_libraries(test_gsl gsl blas)

set(SOURCE_FILES test_gd.cpp gradient_descent.cpp rosenbrock.cpp gradient_descent.h)
add_executable(testGd ${SOURCE_FILES})
target_link_libraries(testGd adept gsl blas)
