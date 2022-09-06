# bb-sbbu
Branch and Bound Algorithm to the Edge Covering Problem on SBBU order.

# Run All Tests

Run all test in parallel using 10 threads. If you like, you could change the number of threads to the one more appropriated to your computational environment.
```
parallel -j 10 < callBB.sh
```

# Development tools

- Intel oneAPI Base Toolkit

# Gtest

my_project$ cmake -S . -B build
-- The C compiler identification is GNU 10.2.1
-- The CXX compiler identification is GNU 10.2.1
...
-- Build files have been written to: .../my_project/build

my_project$ cmake --build build
Scanning dependencies of target gtest
...
[100%] Built target gmock_main

my_project$ cd build && ctest
Test project .../my_project/build
    Start 1: HelloTest.BasicAssertions
1/1 Test #1: HelloTest.BasicAssertions ........   Passed    0.00 sec

100% tests passed, 0 tests failed out of 1

Total Test time (real) =   0.01 sec
Debugging command line
```
cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_C_FLAGS_DEBUG="-g -O0" -DCMAKE_CXX_FLAGS_DEBUG="-g -O0"
```

# References
1. https://blog.ronin.cloud/gnu-parallel/
2. https://google.github.io/googletest/quickstart-cmake.html