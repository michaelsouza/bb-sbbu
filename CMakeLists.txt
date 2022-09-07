cmake_minimum_required(VERSION 3.14)
project(bb-sbbu)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
    tests.bin
    codes/tests.cpp
)

target_link_libraries(
    tests.bin
    GTest::gtest_main
)

add_executable(
  bb.bin
  codes/bb.cpp
)

include(GoogleTest)
gtest_discover_tests(tests.bin)