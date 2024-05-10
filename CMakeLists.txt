cmake_minimum_required(VERSION 3.0.0)
project(cv_takehome_challenge VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

add_executable(cv_takehome_challenge main.cpp)
add_subdirectory(ext/eigen)
add_subdirectory(ext/nlohmann_json)
target_link_libraries(cv_takehome_challenge Eigen3::Eigen nlohmann_json::nlohmann_json)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
