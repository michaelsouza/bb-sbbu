# bb-sbbu
Branch and Bound Algorithm to the Edge Covering Problem on SBBU order.

# Run All Tests

Run all test in parallel using 10 threads. If you like, you could change the number of threads to the one more appropriated to your computational environment.
```
parallel -j 10 < callBB.sh
```

# Development tools

- Visual Studio Code (IDE)
- GCC, the GNU Compiler Collection (compiler)
- GoogleTest (xUnit test framework)
- CMake (Build Framework)
- Git (control version)

# Build Release version

Configure the build in Release mode.
```
> mkdir build

> cmake .. -DCMAKE_BUILD_TYPE=RELEASE

> make
```

# Build Debug version and run the unit tests

Configure the build in Debug mode with all optimizations turned off.
```
> mkdir build

> cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_C_FLAGS_DEBUG="-g -O0" -DCMAKE_CXX_FLAGS_DEBUG="-g -O0"
```
Run the unit tests.
```
> make; ctest
```

# References
1. https://blog.ronin.cloud/gnu-parallel/
2. https://google.github.io/googletest/quickstart-cmake.html