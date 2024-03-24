#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

#include <space/format/Mode.h>

using namespace space;

using D = int;
using V = int;

template <class Mode> void insertRange(Mode &mode, D start, D end) {
  for (D i = start; i < end; i++) {
    mode.insert(i, i);
  }
}

template <class D>
std::vector<D> generateUniform(size_t n, D start, D end, size_t seed = 0) {
  std::vector<D> result;
  std::unordered_set<D> seen;
  std::mt19937 gen(seed);
  std::uniform_int_distribution<D> dis(start, end);
  while (result.size() < n) {
    D x = dis(gen);
    if (seen.count(x) == 0) {
      result.push_back(x);
      seen.insert(x);
    }
  }
  return result;
}

template <class Mode> size_t lookupRange(Mode &mode, D start, D end) {
  size_t count = 0;
  for (D i = start; i < end; i++) {
    if (mode.lookup(i) != nullptr) {
      count++;
    }
  }
  return count;
}

template <class Mode>
size_t lookupVector(Mode &mode, const std::vector<D> &vec) {
  size_t count = 0;
  for (D i : vec) {
    if (mode.lookup(i) != nullptr) {
      count++;
    }
  }
  return count;
}

template <class Mode> size_t iterate(Mode &mode) {
  auto it = mode.iterator();
  size_t count = 0;
  while (it.hasNext()) {
    ++it;
    ++count;
  }
  return count;
}

struct Timer {
  std::chrono::time_point<std::chrono::high_resolution_clock> start;

  Timer() : start(std::chrono::high_resolution_clock::now()) {}

  double tick() {
    auto end = std::chrono::high_resolution_clock::now();
    double result = std::chrono::duration<double>(end - start).count();
    start = end;
    return result;
  }
};

int main() {
  constexpr size_t U = 10000000;
  constexpr double S = 0.3;
  constexpr size_t N = U * S;
  constexpr size_t N2 = N * (S * 0.00001);

  Timer timer;

  std::cout << "== CREATE ==" << std::endl;
  DenseMode<D, D, U> dense_mode;
  std::cout << "DenseMode create: " << timer.tick() << std::endl;
  BitsetDenseMode<D, D, U> bitset_mode;
  std::cout << "BitsetDenseMode create: " << timer.tick() << std::endl;
  SortedMode<D, D, U> sorted_mode;
  std::cout << "SortedMode create: " << timer.tick() << std::endl;
  HashMode<D, D, U> hash_mode;
  std::cout << "HashMode create: " << timer.tick() << std::endl;

  // Insertion - linear
  std::cout << "== INSERT ==" << std::endl;
  insertRange(dense_mode, 0, N);
  std::cout << "DenseMode insert: " << timer.tick() << std::endl;
  insertRange(bitset_mode, 0, N);
  std::cout << "BitsetDenseMode insert: " << timer.tick() << std::endl;
  insertRange(sorted_mode, 0, N);
  std::cout << "SortedMode insert: " << timer.tick() << std::endl;
  insertRange(hash_mode, 0, N);
  std::cout << "HashMode insert: " << timer.tick() << std::endl;

  // Finalize
  std::cout << "== FINALIZE ==" << std::endl;
  sorted_mode.finalize();
  std::cout << "SortedMode finalize: " << timer.tick() << std::endl;

  // Lookup - linear
  size_t r;
  std::cout << "== LOOKUP ==" << std::endl;
  r = lookupRange(dense_mode, 0, N);
  std::cout << "DenseMode lookup: " << r << " in " << timer.tick() << std::endl;
  r = lookupRange(bitset_mode, 0, N);
  std::cout << "BitsetDenseMode lookup: " << r << " in " << timer.tick()
            << std::endl;
  r = lookupRange(sorted_mode, 0, N);
  std::cout << "SortedMode lookup: " << r << " in " << timer.tick()
            << std::endl;
  r = lookupRange(hash_mode, 0, N);
  std::cout << "HashMode lookup: " << r << " in " << timer.tick() << std::endl;

  // Lookup - random
  std::cout << "== LOOKUP RANDOM ==" << std::endl;
  auto vec = generateUniform(N, (D)0, (D)U);
  r = lookupVector(dense_mode, vec);
  std::cout << "DenseMode lookup: " << r << " in " << timer.tick() << std::endl;
  r = lookupVector(bitset_mode, vec);
  std::cout << "BitsetDenseMode lookup: " << r << " in " << timer.tick()
            << std::endl;
  r = lookupVector(sorted_mode, vec);
  std::cout << "SortedMode lookup: " << r << " in " << timer.tick()
            << std::endl;
  r = lookupVector(hash_mode, vec);
  std::cout << "HashMode lookup: " << r << " in " << timer.tick() << std::endl;

  // Iterate
  std::cout << "== ITERATE ==" << std::endl;
  r = iterate(dense_mode);
  std::cout << "DenseMode iterate: " << r << " in " << timer.tick()
            << std::endl;
  r = iterate(bitset_mode);
  std::cout << "BitsetDenseMode iterate: " << r << " in " << timer.tick()
            << std::endl;
  r = iterate(sorted_mode);
  std::cout << "SortedMode iterate: " << r << " in " << timer.tick()
            << std::endl;
  r = iterate(hash_mode);
  std::cout << "HashMode iterate: " << r << " in " << timer.tick() << std::endl;

  // CSR - dense of sorted
  using CSR = BitsetDenseMode<D, SortedMode<D, V, Unbounded>, U>;
  CSR csr;

  // Randomly fill
  auto rows = generateUniform(N, (D)0, (D)U);
  auto cols = generateUniform(N2, (D)0, (D)U);
  std::cout << "Filling with " << N * N2 << " elements" << std::endl;
  timer.tick();
  r = 0;
  for (D row : rows) {
    auto p = csr.insert(row, {});
    for (D col : cols) {
      p->insert(col, 1);
      ++r;
    }
    p->finalize();
  }
  std::cout << "CSR insert: " << timer.tick() << std::endl;

  return 0;
}