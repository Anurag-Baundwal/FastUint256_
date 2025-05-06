#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include <random> 

// For Boost
#include <boost/multiprecision/cpp_int.hpp>

// --- Include FastUint256 header ---
#include "fast_uint256.h" // Assuming it's in the same directory or include path

// --- We use the namespace defined in the FastUint256 header ---
using namespace my_math;

// Define the Boost type alias
using BoostUint256 = boost::multiprecision::uint256_t;

// --- Utility Functions for Random Data ---
// Generate a random 256-bit hex string (e.g., "0x...")
std::string generate_random_hex256(std::mt19937_64 &rng) {
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  std::stringstream ss;
  ss << "0x";
  bool all_zero = true;
  std::array<uint64_t, 4> limbs_val;
  for (int i = 0; i < 4; ++i) {
    limbs_val[i] = dist(rng);
    if (limbs_val[i] != 0) all_zero = false;
  }
  if (all_zero) return "0x0"; // Handle zero case explicitly
  bool leading_zeros = true;
  for (int i = 3; i >= 0; --i) {
    if (limbs_val[i] == 0 && leading_zeros && i > 0) continue;
    if (leading_zeros) {
      ss << std::hex << limbs_val[i];
      leading_zeros = false;
    } else {
      ss << std::hex << std::setfill('0') << std::setw(16) << limbs_val[i];
    }
  }
  return ss.str();
}

// Convert BoostUint256 to a hex string
// (This is a helper for potentially printing Boost results, not directly used
// by FastUint256)
std::string boost_to_hex_string(const BoostUint256 &val,
                                bool include_prefix = true) {
  if (val == 0) return include_prefix ? "0x0" : "0";
  std::stringstream ss;
  if (include_prefix) ss << "0x";
  ss << std::hex << val;
  return ss.str();
}

// --- Main Benchmarking Function ---
int main() {
  const size_t NUM_PAIRS = 10; // Number of unique data pairs to cycle through
  const size_t NUM_REPETITIONS =
      100000; // Number of times to repeat operations on the set of pairs
  const size_t TOTAL_OPS =
      NUM_PAIRS * NUM_REPETITIONS; // Total operations per benchmark section

  // Setup random number generation
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<uint64_t> scalar_dist(
      0, std::numeric_limits<uint64_t>::max());
  std::uniform_int_distribution<size_t> shift_dist(0, 255); // Shifts from 0 to 255

  std::cout << "--- Starting Benchmark ---" << std::endl;
  std::cout << "Operations per benchmark section: " << TOTAL_OPS << std::endl;

  // --- Print info based on macros defined IN fast_uint256.h ---
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
  std::cout << "----------------------------------------" << std::endl;

  // --- Generate Test Data ---
  std::vector<std::pair<FastUint256, FastUint256>> data_fu(NUM_PAIRS);
  std::vector<std::pair<BoostUint256, BoostUint256>> data_bu(NUM_PAIRS);
  std::vector<uint64_t> scalars(NUM_PAIRS);
  std::vector<size_t> shifts(NUM_PAIRS);

  std::cout << "Generating random test data..." << std::endl;
  for (size_t i = 0; i < NUM_PAIRS; ++i) {
    std::string hex_a = generate_random_hex256(rng);
    std::string hex_b = generate_random_hex256(rng);

    try {
      data_fu[i].first = FastUint256(hex_a);
      data_fu[i].second = FastUint256(hex_b);

      std::stringstream ss_a(hex_a);
      ss_a >> std::hex >> data_bu[i].first;
      if (ss_a.fail() && !ss_a.eof())
        throw std::runtime_error("Boost parse A failed");
      std::stringstream ss_b(hex_b);
      ss_b >> std::hex >> data_bu[i].second;
      if (ss_b.fail() && !ss_b.eof())
        throw std::runtime_error("Boost parse B failed");

      scalars[i] = scalar_dist(rng);
      shifts[i] = shift_dist(rng);
    } catch (const std::exception &e) {
      std::cerr << "ERROR generating data for index " << i << ": " << e.what()
                << std::endl;
      return 1;
    }
  }
  std::cout << "Data generation complete." << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmarking Variables ---
  std::chrono::high_resolution_clock::time_point start_time, end_time;
  double duration_fu_ms, duration_bu_ms;

  // --- Benchmark Bitwise AND (&) ---
  std::cout << "Benchmarking: Bitwise AND (&)" << std::endl;
  FastUint256 result_fu_and = {}; // Accumulator
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_and ^= data_fu[j].first & data_fu[j].second; // Accumulate
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_and = 0; // Accumulator
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_and ^= data_bu[j].first & data_bu[j].second; // Accumulate
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  // Optional: Print sink to ensure usage (and correctness if values were known)
  // std::cout << "  (Sink FU: " << result_fu_and << ", Sink BU: " <<
  // boost_to_hex_string(result_bu_and) << ")" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Bitwise OR (|) ---
  std::cout << "Benchmarking: Bitwise OR (|)" << std::endl;
  FastUint256 result_fu_or = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_or ^= data_fu[j].first | data_fu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_or = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_or ^= data_bu[j].first | data_bu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Bitwise XOR (^) ---
  std::cout << "Benchmarking: Bitwise XOR (^)" << std::endl;
  FastUint256 result_fu_xor = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_xor ^= data_fu[j].first ^ data_fu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_xor = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_xor ^= data_bu[j].first ^ data_bu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Bitwise NOT (~) ---
  std::cout << "Benchmarking: Bitwise NOT (~)" << std::endl;
  FastUint256 result_fu_not = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_not ^= ~data_fu[j].first; // Only use first element of pair
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_not = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_not ^= ~data_bu[j].first; // Only use first element of pair
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Left Shift (<<) ---
  std::cout << "Benchmarking: Left Shift (<<)" << std::endl;
  FastUint256 result_fu_shl = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_shl ^= data_fu[j].first << shifts[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_shl = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_shl ^= data_bu[j].first << shifts[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Right Shift (>>) ---
  std::cout << "Benchmarking: Right Shift (>>)" << std::endl;
  FastUint256 result_fu_shr = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_shr ^= data_fu[j].first >> shifts[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_shr = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_shr ^= data_bu[j].first >> shifts[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Addition (+) ---
  std::cout << "Benchmarking: Addition (+)" << std::endl;
  FastUint256 result_fu_add = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_add ^= data_fu[j].first + data_fu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_add = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_add ^= data_bu[j].first + data_bu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Subtraction (-) ---
  std::cout << "Benchmarking: Subtraction (-)" << std::endl;
  FastUint256 result_fu_sub = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_sub ^= data_fu[j].first - data_fu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_sub = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_sub ^= data_bu[j].first - data_bu[j].second;
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  // --- Benchmark Scalar Multiplication (*) ---
  std::cout << "Benchmarking: Scalar Multiplication (* uint64_t)" << std::endl;
  FastUint256 result_fu_mul = {};
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_fu_mul ^= data_fu[j].first * scalars[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_fu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  FastUint256: " << std::fixed << std::setprecision(2)
            << duration_fu_ms << " ms" << std::endl;

  BoostUint256 result_bu_mul = 0;
  start_time = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < NUM_REPETITIONS; ++i) {
    for (size_t j = 0; j < NUM_PAIRS; ++j) {
      result_bu_mul ^= data_bu[j].first * scalars[j];
    }
  }
  end_time = std::chrono::high_resolution_clock::now();
  duration_bu_ms =
      std::chrono::duration<double, std::milli>(end_time - start_time).count();
  std::cout << "  BoostUint256:  " << std::fixed << std::setprecision(2)
            << duration_bu_ms << " ms" << std::endl;
  std::cout << "----------------------------------------" << std::endl;

  std::cout << "--- Benchmark Complete ---" << std::endl;

  return 0;
}