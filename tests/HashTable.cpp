#include <space/format/HashTable.h>

#include <array>
#include <chrono>
#include <iostream>
#include <random>
#include <space/util/Fn.h>
#include <tuple>

// WARNING: not a benchmark! just making sure things kinda work

using Key = uint64_t;
using Value = int64_t;

using namespace space;

bool test_1() {
  constexpr size_t N = 10000000;

  HashTable<Key, Value, /*Instrument=*/true> ht(1, 1, 2);

  // randomly generate keys and values
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Key> dis(0, N);

  std::vector<std::tuple<Key>> keys;

  // time insert
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    Key k[1] = {i};
    keys.push_back(std::make_tuple(k[0]));
    ht.insert(k, i);
  }
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();

  double insert_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  // time lookup
  t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    auto [k0] = keys[i];
    Key k[1] = {k0};
    Value *v = ht.lookup(k);
    if (!v) {
      std::cerr << "key not found: " << k0 << std::endl;
      return false;
    }
    if (static_cast<size_t>(*v) != i) {
      std::cerr << "wrong value: " << k0 << std::endl;
      std::cerr << "  expected: " << i << std::endl;
      std::cerr << "  got: " << *v << std::endl;
      return false;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();

  double lookup_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  std::cout << "HashTable test 1 passed" << std::endl;
  std::cout << "  lookups: " << ht.lookups << std::endl;
  std::cout << "  probes: " << ht.probes << std::endl;
  std::cout << "  probes / lookup: " << (double)ht.probes / ht.lookups
            << std::endl;
  std::cout << "  final capacity: " << ht.getCapacity() << " ("
            << (double)ht.getSize() / ht.getCapacity() * 100 << "% full)"
            << std::endl;
  std::cout << "  memory usage: " << ht.getMemoryUsage() / (1024 * 1024)
            << " MB" << std::endl;
  std::cout << "  insert time: " << insert_time_ms << " ms" << std::endl;
  std::cout << "  lookup time: " << lookup_time_ms << " ms" << std::endl;

  return true;
}

bool test_2() {
  constexpr size_t N = 10000000;

  // randomly generate keys and values
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<Key> dis(0, N);

  std::vector<std::tuple<Key, Key>> keys;

  // time insert
  HashTable<Key, Value, /*Instrument=*/true> ht(2, 1, 3);
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    Key k[2] = {i, i};
    // while (ht.lookup(k)) {
    //   k[0] = dis(gen);
    //   k[1] = dis(gen);
    // }
    keys.push_back(std::make_tuple(k[0], k[1]));
    ht.insert(k, i);
  }
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();

  double insert_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (ht.getSize() != N) {
    std::cerr << "wrong size: " << ht.getSize() << std::endl;
    return false;
  }

  // time lookup
  t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    auto [k0, k1] = keys[i];
    Key k[2] = {k0, k1};
    Value *v = ht.lookup(k);
    if (!v) {
      std::cerr << "key not found: " << k0 << ", " << k1 << std::endl;
      return false;
    }
    if (static_cast<size_t>(*v) != i) {
      std::cerr << "wrong value: " << k0 << ", " << k1 << std::endl;
      return false;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();

  double lookup_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  // time iter
  std::array<Key, 2> out;
  Value v;
  t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    ht.iter(i, out.data(), &v);
    if (v != static_cast<Value>(out[0])) {
      std::cerr << "wrong value: " << out[0] << std::endl;
      return false;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();

  double iter_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  std::cout << "HashTable test 2 passed" << std::endl;
  std::cout << "  lookups: " << ht.lookups << std::endl;
  std::cout << "  probes: " << ht.probes << std::endl;
  std::cout << "  probes / lookup: " << (double)ht.probes / ht.lookups
            << std::endl;
  std::cout << "  final capacity: " << ht.getCapacity() << " ("
            << (double)ht.getSize() / ht.getCapacity() * 100 << "% full)"
            << std::endl;
  std::cout << "  memory usage: " << ht.getMemoryUsage() / (1024 * 1024)
            << " MB" << std::endl;
  std::cout << "  insert time: " << insert_time_ms << " ms" << std::endl;
  std::cout << "  lookup time: " << lookup_time_ms << " ms" << std::endl;
  std::cout << "  iter time: " << iter_time_ms << " ms" << std::endl;

  return true;
}

bool test_2_static() {
  constexpr size_t N = 10000000;

  // std::vector<std::tuple<Key, Key>> keys;
  std::vector<std::array<Key, 2>> keys;

  // time insert
  // HashTable<Key, Value, /*Instrument=*/true> ht(2, 1, 3);
  StaticHashTable<Key, Value, 2> ht(1, 3);
  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    keys.push_back({i, i});
    ht.insert(keys[i], i);
  }
  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();

  double insert_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (ht.getSize() != N) {
    std::cerr << "wrong size: " << ht.getSize() << std::endl;
    return false;
  }

  // time lookup
  t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    auto [k0, k1] = keys[i];
    Value *v = ht.lookup(keys[i]);
    if (!v) {
      std::cerr << "key not found: " << k0 << ", " << k1 << std::endl;
      return false;
    }
    if (static_cast<size_t>(*v) != i) {
      std::cerr << "wrong value: " << k0 << ", " << k1 << std::endl;
      return false;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();

  double lookup_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  // time iter
  std::array<Key, 2> out;
  Value v;
  t1 = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < N; i++) {
    ht.iter(i, out, v);
    if (v != static_cast<Value>(out[0])) {
      std::cerr << "wrong value: " << out[0] << std::endl;
      return false;
    }
  }
  t2 = std::chrono::high_resolution_clock::now();

  double iter_time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  std::cout << "HashTable test 2 passed" << std::endl;
  // std::cout << "  lookups: " << ht.lookups << std::endl;
  // std::cout << "  probes: " << ht.probes << std::endl;
  // std::cout << "  probes / lookup: " << (double)ht.probes / ht.lookups
  //           << std::endl;
  std::cout << "  final capacity: " << ht.getCapacity() << " ("
            << (double)ht.getSize() / ht.getCapacity() * 100 << "% full)"
            << std::endl;
  std::cout << "  memory usage: " << ht.getMemoryUsage() / (1024 * 1024)
            << " MB" << std::endl;
  std::cout << "  insert time: " << insert_time_ms << " ms" << std::endl;
  std::cout << "  lookup time: " << lookup_time_ms << " ms" << std::endl;
  std::cout << "  iter time: " << iter_time_ms << " ms" << std::endl;

  return true;
}

int main() {
  if (!test_1()) {
    return 1;
  }
  if (!test_2()) {
    return 1;
  }
  if (!test_2_static()) {
    return 1;
  }

  return 0;
}
