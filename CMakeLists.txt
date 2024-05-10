cmake_minimum_required(VERSION 3.0.0)
project(cv_takehome_challenge VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(cv_takehome_challenge main.cpp)
find_package(Eigen3 3.4 NO_MODULE)
target_link_libraries(cv_takehome_challenge Eigen3::Eigen)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
