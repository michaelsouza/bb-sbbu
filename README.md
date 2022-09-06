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

Build in Debug mode with all optimization turned off
```
cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_C_FLAGS_DEBUG="-g -O0" -DCMAKE_CXX_FLAGS_DEBUG="-g -O0"
```

# References
1. https://blog.ronin.cloud/gnu-parallel/
2. https://google.github.io/googletest/quickstart-cmake.html