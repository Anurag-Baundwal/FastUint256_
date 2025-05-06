# FastUint256

A C++ header-only library providing a `FastUint256` type for 256-bit unsigned integer operations.

## Files

*   `FastUint256.h`: The header library.
*   `FastUint256BoostTest.cpp`: Correctness tests against `boost::multiprecision::uint256_t`.
*   `FastUint256Benchmark.cpp`: Performance benchmarks against `boost::multiprecision::uint256_t`.

## Usage

1.  Include `FastUint256.h`.
2.  Use the type `my_math::FastUint256`.

```cpp
#include "FastUint256.h" // Or <path/to/FastUint256.h>
#include <iostream>

int main() {
    my_math::FastUint256 a("0x1234567890abcdef1234567890abcdef");
    my_math::FastUint256 b(1000);
    std::cout << "A + B = " << (a + b) << std::endl;
    return 0;
}
```

## Building Tests & Benchmarks

The test and benchmark executables require the Boost C++ Libraries to use `boost::multiprecision::uint256_t` as a reference for correctness and performance. The `FastUint256.h` library itself does not depend on Boost.

Example (MSVC):

### Test
`cl /std:c++17 /EHsc /O2 /I C:\path\to\boost\include FastUint256BoostTest.cpp /Fe:FastUint256BoostTest.exe`

### Benchmark
`cl /std:c++17 /EHsc /O2 /I C:\path\to\boost\include FastUint256Benchmark.cpp /Fe:FastUint256Benchmark.exe /DNDEBUG`

(Replace `C:\path\to\boost\include` with your Boost path (on Windows this is usually `C:\Boost\boost_1_84_0` or similar, depending on the version number))