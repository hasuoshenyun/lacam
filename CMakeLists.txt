cmake_minimum_required(VERSION 3.16)
project(lacam-project CXX)

# 全局C++标准设置
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 设置编译器路径，仅在未设置时更新缓存
#set(CMAKE_C_COMPILER /usr/bin/clang CACHE STRING "C compiler" FORCE)
#set(CMAKE_CXX_COMPILER /usr/bin/clang++ CACHE STRING "C++ compiler" FORCE)

set(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang CACHE STRING "C compiler" FORCE)
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++ CACHE STRING "C++ compiler" FORCE)
include_directories(/opt/homebrew/opt/llvm/include/c++/v1)

# 设置调试模式
set(CMAKE_BUILD_TYPE Debug)  # 强制使用调试模式

# 禁用优化（以便调试时能单步调试）
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")  # -g 用于生成调试信息，-O0 禁用优化

# 打印调试信息
message(STATUS "C++ Standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

# 添加子模块
add_subdirectory(./lacam)
add_subdirectory(./third_party/argparse)
add_subdirectory(./third_party/googletest)


# 主程序
add_executable(main main.cpp)
target_compile_features(main PUBLIC cxx_std_17)
target_link_libraries(main lacam argparse)

# 测试设置
set(TEST_MAIN_FUNC ./third_party/googletest/googletest/src/gtest_main.cc)
set(TEST_ALL_SRC ${TEST_MAIN_FUNC})

# 添加测试宏
macro(add_test name target)
  add_executable(${name} ${target} ${TEST_MAIN_FUNC})
  target_link_libraries(${name} lacam gtest)
  list(APPEND TEST_ALL_SRC ${target})
endmacro(add_test)

# 定义测试用例
add_test(test_graph ./tests/test_graph.cpp)
add_test(test_instance ./tests/test_instance.cpp)
add_test(test_dist_table ./tests/test_dist_table.cpp)
add_test(test_planner ./tests/test_planner.cpp)
add_test(test_post_processing ./tests/test_post_processing.cpp)

# 汇总所有测试
add_executable(test_all ${TEST_ALL_SRC})
target_link_libraries(test_all lacam gtest)

