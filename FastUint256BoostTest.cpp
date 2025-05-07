// Core Standard Libraries needed for Testing
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include <limits>    // For numeric_limits
#include <stdexcept> // For exception handling during tests
#include <cmath>     // For std::abs in potential future float tests (not needed here)
#include <cassert>   // For potential assertions (alternative to cout checks)

// Boost Multiprecision
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/detail/bitscan.hpp> // For msb/lsb if needed directly

// --- Include FastUint256 header ---
#include "FastUint256.h"

// --- We use the namespace defined in the header ---
using namespace my_math;

// ========================================
// --- Start of Testing Infrastructure ---
// ========================================

// Define the Boost type alias
using BoostUint256 = boost::multiprecision::uint256_t;

// --- Utility Functions for Testing ---

// Generate a random 256-bit hex string (e.g., "0x...")
std::string generate_random_hex256(std::mt19937_64 &rng) {
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  std::stringstream ss;
  ss << "0x";
  bool all_zero = true;
  std::array<uint64_t, 4> limbs;
  for (int i = 0; i < 4; ++i) {
    limbs[i] = dist(rng);
    if (limbs[i] != 0) all_zero = false;
  }

  if (all_zero) return "0x0"; // Handle zero case explicitly

  bool leading_zeros = true;
  for (int i = 3; i >= 0; --i) {
    if (limbs[i] == 0 && leading_zeros && i > 0) {
      continue; // Skip leading zero limbs unless it's the only one
    }
    if (leading_zeros) {
      // Print first non-zero limb without padding
      ss << std::hex << limbs[i];
      leading_zeros = false;
    } else {
      // Print subsequent limbs with full 16-char padding
      ss << std::hex << std::setfill('0') << std::setw(16) << limbs[i];
    }
  }
  return ss.str();
}

// Convert BoostUint256 to a hex string compatible with
// FastUint256::to_hex_string
std::string boost_to_hex_string(const BoostUint256 &val,
                                bool include_prefix = true) {
  if (val == 0) {
    return include_prefix ? "0x0" : "0";
  }
  std::stringstream ss;
  if (include_prefix) {
    ss << "0x";
  }
  ss << std::hex << val;
  return ss.str();
}

// Comparison function for results (FastUint256 vs BoostUint256)
template <typename T_Fast, typename T_Boost>
bool check_results(const std::string &op_name, const std::string &hex_a,
                   const std::string &hex_b, // Inputs or context
                   const T_Fast &result_fu, const T_Boost &result_bu,
                   int &fail_count) {
  std::string fu_hex = result_fu.to_hex_string();
  std::string bu_hex = boost_to_hex_string(result_bu);

  if (fu_hex == bu_hex) {
    return true;
  } else {
    std::cout << "  " << op_name << ": FAIL" << std::endl;
    std::cout << "    Operand A: " << hex_a << std::endl;
    if (!hex_b.empty()) { // Display second operand if provided (e.g., B,
                          // scalar, shift)
      std::cout << "    Operand B/Param: " << hex_b << std::endl;
    }
    std::cout << "    FastUint256 Result: " << fu_hex << std::endl;
    std::cout << "    BoostUint256 Result: " << bu_hex << std::endl;
    fail_count++;
    return false;
  }
}

// Comparison function for boolean results
bool check_bool_result(const std::string &op_name, const std::string &hex_a,
                       const std::string &hex_b, // Inputs
                       bool result_fu, bool result_bu, int &fail_count) {
  if (result_fu == result_bu) {
    return true;
  } else {
    std::cout << "  " << op_name << ": FAIL" << std::endl;
    std::cout << "    Input A: " << hex_a << std::endl;
    std::cout << "    Input B: " << hex_b << std::endl;
    std::cout << "    FastUint256: " << (result_fu ? "true" : "false")
              << std::endl;
    std::cout << "    BoostUint256: " << (result_bu ? "true" : "false")
              << std::endl;
    fail_count++;
    return false;
  }
}

// Comparison function for integer results (clz, ctz, popcount)
bool check_int_result(const std::string &op_name,
                      const std::string &hex_val, // Input
                      int result_fu, int result_bu, int &fail_count) {
  if (result_fu == result_bu) {
    return true;
  } else {
    std::cout << "  " << op_name << ": FAIL" << std::endl;
    std::cout << "    Input: " << hex_val << std::endl;
    std::cout << "    FastUint256: " << result_fu << std::endl;
    std::cout << "    BoostUint256: " << result_bu << std::endl;
    fail_count++;
    return false;
  }
}

// Calculate clz for BoostUint256
int boost_clz(const BoostUint256 &val) {
  if (val == 0) {
    return 256;
  }
  unsigned msb_index = boost::multiprecision::msb(val);
  return 255 - msb_index;
}

// Calculate ctz for BoostUint256
int boost_ctz(const BoostUint256 &val) {
  if (val == 0) {
    return 256;
  }
  return boost::multiprecision::lsb(val);
}

// Calculate popcount for BoostUint256
int boost_popcount(const BoostUint256 &val) {
  if (val == 0) {
    return 0;
  }
  // Use Kernighan's algorithm (efficient for sparse bits)
  // BoostUint256 supports the necessary bitwise and arithmetic ops
  BoostUint256 temp = val; // Work on a copy
  int count = 0;
  while (temp != 0) {
    temp &= (temp - 1); // Clear the least significant bit set
    count++;
  }
  return count;

  // Alternative
  // Iterate bits
  // int count = 0;
  // for(unsigned i = 0; i < 256; ++i) {
  //     if (boost::multiprecision::bit_test(val, i)) {
  //         count++;
  //     }
  // }
  // return count;
}

// --- Main Testing Function ---
int main() {
  const int NUM_TESTS = 100; // Number of random pairs to test
  int total_failures = 0;

  // Setup random number generation
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<uint64_t> scalar_dist(
      0, std::numeric_limits<uint64_t>::max());
  std::uniform_int_distribution<size_t> shift_dist(0,
                                                   256); // Include shift by 256

  std::cout << "--- Starting Boost Comparison Test (" << NUM_TESTS
            << " random pairs) ---" << std::endl;

  // Print info based on macros defined IN fast_uint256.h
#if FASTUINT_HAS_MSVC_INTRINSICS && defined(_M_X64)
  std::cout << "INFO: FastUint256 using MSVC Intrinsics." << std::endl;
#elif FASTUINT_HAS_GCC_CLANG_INTRINSICS && defined(__x86_64__) &&              \
    (defined(__ADX__) || defined(__POPCNT__) || defined(__LZCNT__) ||          \
     defined(__BMI__))
  std::cout << "INFO: FastUint256 using GCC/Clang Intrinsics "
               "(ADX/POPCNT/LZCNT/BMI detected)."
            << std::endl;
#elif FASTUINT_HAS_UINT128
  std::cout << "INFO: FastUint256 using __uint128_t." << std::endl;
#else
  std::cout << "INFO: FastUint256 using Portable Implementation." << std::endl;
#endif

  for (int i = 0; i < NUM_TESTS; ++i) {
    std::string hex_a = generate_random_hex256(rng);
    std::string hex_b = generate_random_hex256(rng);
    uint64_t scalar = scalar_dist(rng);
    size_t shift = shift_dist(rng);
    // Add some edge case shift values occasionally
    if (i % 10 == 0) shift = 0;
    if (i % 10 == 1) shift = 63;
    if (i % 10 == 2) shift = 64;
    if (i % 10 == 3) shift = 128;
    if (i % 10 == 4) shift = 255;
    if (i % 10 == 5) shift = 256;

    int test_case_failures = 0;
    std::cout << "\n--- Test Case " << (i + 1) << " ---" << std::endl;
    // Print test case details
    // std::cout << "  A = " << hex_a << std::endl;
    // std::cout << "  B = " << hex_b << std::endl;
    // std::cout << "  Scalar = " << scalar << std::endl;
    // std::cout << "  Shift = " << shift << std::endl;

    try {
      // --- Create instances ---
      FastUint256 a_fu(hex_a);
      FastUint256 b_fu(hex_b);
      BoostUint256 a_bu;
      BoostUint256 b_bu;

      // Boost construction
      try {
        std::stringstream ss_a(hex_a);
        ss_a >> std::hex >> a_bu;
        if (ss_a.fail() && !ss_a.eof())
          throw std::runtime_error("Boost failed to parse hex A");
        std::stringstream ss_b(hex_b);
        ss_b >> std::hex >> b_bu;
        if (ss_b.fail() && !ss_b.eof())
          throw std::runtime_error("Boost failed to parse hex B");
      } catch (const std::exception &e) {
        std::cerr << "ERROR: Boost construction failed for A=" << hex_a
                  << " or B=" << hex_b << ": " << e.what() << std::endl;
        total_failures++;
        continue;
      }

      // --- Perform Operations & Compare ---

      // --- Bitwise ---
      check_results("AND (&)", hex_a, hex_b, a_fu & b_fu, a_bu & b_bu,
                    test_case_failures);
      { // Test &=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu &= b_fu;
        temp_a_bu &= b_bu;
        check_results("AND Assign (A &= B)", hex_a, hex_b, temp_a_fu, temp_a_bu,
                      test_case_failures);
      }

      check_results("OR (|)", hex_a, hex_b, a_fu | b_fu, a_bu | b_bu,
                    test_case_failures);
      { // Test |=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu |= b_fu;
        temp_a_bu |= b_bu;
        check_results("OR Assign (A |= B)", hex_a, hex_b, temp_a_fu, temp_a_bu,
                      test_case_failures);
      }

      check_results("XOR (^)", hex_a, hex_b, a_fu ^ b_fu, a_bu ^ b_bu,
                    test_case_failures);
      { // Test ^=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu ^= b_fu;
        temp_a_bu ^= b_bu;
        check_results("XOR Assign (A ^= B)", hex_a, hex_b, temp_a_fu, temp_a_bu,
                      test_case_failures);
      }

      check_results("NOT (~A)", hex_a, "", ~a_fu, ~a_bu, test_case_failures);
      // No compound assignment for ~

      // --- Shifts ---
      std::string shift_str = std::to_string(shift);
      check_results("Left Shift (A << S)", hex_a, shift_str, a_fu << shift,
                    a_bu << shift, test_case_failures);
      { // Test <<=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu <<= shift;
        temp_a_bu <<= shift;
        check_results("Left Shift Assign (A <<= S)", hex_a, shift_str,
                      temp_a_fu, temp_a_bu, test_case_failures);
      }

      check_results("Right Shift (A >> S)", hex_a, shift_str, a_fu >> shift,
                    a_bu >> shift, test_case_failures);
      { // Test >>=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu >>= shift;
        temp_a_bu >>= shift;
        check_results("Right Shift Assign (A >>= S)", hex_a, shift_str,
                      temp_a_fu, temp_a_bu, test_case_failures);
      }

      // --- Arithmetic ---
      check_results("Addition (A + B)", hex_a, hex_b, a_fu + b_fu, a_bu + b_bu,
                    test_case_failures);
      { // Test +=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu += b_fu;
        temp_a_bu += b_bu;
        check_results("Add Assign (A += B)", hex_a, hex_b, temp_a_fu, temp_a_bu,
                      test_case_failures);
      }

      check_results("Subtraction (A - B)", hex_a, hex_b, a_fu - b_fu,
                    a_bu - b_bu, test_case_failures);
      { // Test -=
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu -= b_fu;
        temp_a_bu -= b_bu;
        check_results("Subtract Assign (A -= B)", hex_a, hex_b, temp_a_fu,
                      temp_a_bu, test_case_failures);
      }

      std::string scalar_str = std::to_string(scalar);
      check_results("Scalar Mul (A * scalar)", hex_a, scalar_str, a_fu * scalar,
                    a_bu * scalar, test_case_failures);
      { // Test *= scalar
        FastUint256 temp_a_fu = a_fu;
        BoostUint256 temp_a_bu = a_bu;
        temp_a_fu *= scalar;
        temp_a_bu *= scalar;
        check_results("Scalar Mul Assign (A *= scalar)", hex_a, scalar_str,
                      temp_a_fu, temp_a_bu, test_case_failures);
      }

      // --- Comparisons ---
      check_bool_result("Equal (A == B)", hex_a, hex_b, a_fu == b_fu,
                        a_bu == b_bu, test_case_failures);
      check_bool_result("Not Equal (A != B)", hex_a, hex_b, a_fu != b_fu,
                        a_bu != b_bu, test_case_failures);
      check_bool_result("Less Than (A < B)", hex_a, hex_b, a_fu < b_fu,
                        a_bu < b_bu, test_case_failures);
      check_bool_result("Greater Than (A > B)", hex_a, hex_b, a_fu > b_fu,
                        a_bu > b_bu, test_case_failures);
      check_bool_result("Less Equal (A <= B)", hex_a, hex_b, a_fu <= b_fu,
                        a_bu <= b_bu, test_case_failures);
      check_bool_result("Greater Equal (A >= B)", hex_a, hex_b, a_fu >= b_fu,
                        a_bu >= b_bu, test_case_failures);

      // --- Utility Functions ---
      check_int_result("Count Leading Zeros (clz A)", hex_a, a_fu.clz(),
                       boost_clz(a_bu), test_case_failures);
      check_int_result("Count Trailing Zeros (ctz A)", hex_a, a_fu.ctz(),
                       boost_ctz(a_bu), test_case_failures);

      // --- Popcount ---
      check_int_result("Popcount (A)", hex_a, a_fu.popcount(),
                       boost_popcount(a_bu), test_case_failures);
      // Also check popcount for B for more coverage
      check_int_result("Popcount (B)", hex_b, b_fu.popcount(),
                       boost_popcount(b_bu), test_case_failures);

      // --- Test other utilities ---
      // Test is_zero() and operator bool() ---
      check_bool_result("Is Zero (A.is_zero())", hex_a, "", a_fu.is_zero(),
                        a_bu == 0, test_case_failures);
      check_bool_result("Operator bool (A)", hex_a, "", static_cast<bool>(a_fu),
                        static_cast<bool>(a_bu != 0),
                        test_case_failures); // Boost doesn't have direct
                                             // operator bool, compare with != 0

      // Test operator uint64_t() - only compares if value fits in uint64_t
      // Note: Boost conversion throws if value > uint64_t::max, FastUint just
      // truncates. This test only makes sense if the value fits in 64 bits.
      // Let's check explicitly.
      if (a_bu <= std::numeric_limits<uint64_t>::max()) {
        check_int_result("Operator uint64_t (A)", hex_a,
                         static_cast<uint64_t>(a_fu),
                         a_bu.convert_to<uint64_t>(), test_case_failures);
      } // Otherwise, the behaviour differs by design (truncation vs exception)

    } catch (const std::exception &e) {
      std::cerr << "ERROR during Test Case " << (i + 1) << " (A=" << hex_a
                << ", B=" << hex_b << "): " << e.what() << std::endl;
      test_case_failures++; // Count exception as a failure
    }

    if (test_case_failures == 0) {
      std::cout << "--- Test Case " << (i + 1) << ": ALL PASSED ---"
                << std::endl;
    } else {
      std::cout << "--- Test Case " << (i + 1) << ": FAILED ("
                << test_case_failures << " errors) ---" << std::endl;
      total_failures += test_case_failures; // Accumulate failures
    }
  } // End of test loop

  std::cout << "\n--- Boost Comparison Test Summary ---" << std::endl;
  if (total_failures == 0) {
    std::cout << "Result: ALL " << NUM_TESTS << " test cases PASSED."
              << std::endl;
  } else {
    std::cout << "Result: FAILED - " << total_failures
              << " errors detected across all test cases." << std::endl;
  }

  return (total_failures == 0) ? 0 : 1; // Return 0 on success, 1 on failure
}