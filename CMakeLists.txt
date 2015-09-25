cmake_minimum_required(VERSION 3.3)
project(filterSqp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES test_ad.cpp rosenbrock.cpp)
add_executable(testAd ${SOURCE_FILES})
target_link_libraries(testAd)

add_executable(test_gsl test_gsl.cpp)
target_link_libraries(test_gsl gsl blas)

set(GD_SOURCE_FILES test_gd.cpp gradient_descent.cpp rosenbrock.cpp gradient_descent.h)
add_executable(testGd ${GD_SOURCE_FILES})
target_link_libraries(testGd adept gsl blas gc)

set(FB_SOURCE_FILES test_fadbad.cpp)
add_executable(testFadbad ${FB_SOURCE_FILES})
#target_link_libraries(testFadbad)
