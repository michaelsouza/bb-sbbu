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

# Find GMP and GMPXX libraries
find_library(GMP_LIBRARY gmp)
find_library(GMPXX_LIBRARY gmpxx)

target_link_libraries(
    tests.bin
    GTest::gtest_main
    ${GMP_LIBRARY}
    ${GMPXX_LIBRARY}
)

target_link_libraries(
  bb.bin
  ${GMP_LIBRARY}
  ${GMPXX_LIBRARY}
)

include(GoogleTest)
gtest_discover_tests(tests.bin)