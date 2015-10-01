cmake_minimum_required(VERSION 3.3)
project(filterSqp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

set(SOURCE_FILES test_ad.cpp rosenbrock.cpp)
add_executable(testAd ${SOURCE_FILES})
target_link_libraries(testAd)

add_executable(test_gsl test_gsl.cpp)
target_link_libraries(test_gsl gsl blas)

set(GD_SOURCE_FILES test_gd.cpp gradient_descent.cpp rosenbrock.cpp gradient_descent.h tr_optim.h)
add_executable(testGd ${GD_SOURCE_FILES})
target_link_libraries(testGd gsl blas)

set(LM_SOURCE_FILES test_lm.cpp tr_optim.cpp rosenbrock.cpp himmelblau.cpp tr_optim.h find_root.cpp himmelblau.h)
add_executable(testLm ${LM_SOURCE_FILES})
target_link_libraries(testLm gsl blas)

set(FB_SOURCE_FILES test_fadbad.cpp)
add_executable(testFadbad ${FB_SOURCE_FILES})
#target_link_libraries(testFadbad)

set(TR_SOURCE_FILES test_findroot.cpp find_root.cpp)
add_executable(testFindroot ${TR_SOURCE_FILES})

set(LMPATO_SOURCE_FILES test_lm_pato.cpp tr_optim.cpp find_root.cpp)
add_executable(testLmPato ${LMPATO_SOURCE_FILES})
target_link_libraries(testLmPato gsl blas)

set(PAPER_SOURCE_FILES test_paper.cpp papersheet.cpp tr_optim.cpp find_root.cpp)
add_executable(testPaper ${PAPER_SOURCE_FILES})
target_link_libraries(testPaper gsl blas)

set(PAPER_SOURCE_FILES test_snakeA.cpp snakeA.cpp tr_optim.cpp find_root.cpp)
add_executable(testSnakeA ${PAPER_SOURCE_FILES})
target_link_libraries(testSnakeA gsl blas gc)
