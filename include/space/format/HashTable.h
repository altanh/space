#pragma once

#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>

namespace space {

// We want a hash table from (x1, ..., xn) -> V, where each xi is an integer.
// It is imperative for locality that the key tuples be stored packed in memory,
// rather than pointers to keys. This is because the keys are often small.
//
// This is easy to do statically, but we want to support dynamic tuple sizes due
// to runtime planning. Therefore, for a hash table of size N with key size K,
// we store the entries in a single array of length N * K.
//
// Ideally we would JIT this stuff...

// NB: insert only
template <class Key, class Value, bool Instrument = false> class HashTable {
  static constexpr Key EMPTY_KEY = std::numeric_limits<Key>::max();
  static constexpr size_t INITIAL_CAPACITY = 16;
  static constexpr uint8_t DEFAULT_LF_NUM = 1;
  static constexpr uint8_t DEFAULT_LF_DEN = 2;

  uint8_t tuple_size;
  const size_t entry_size;
  uint8_t lf_num;
  uint8_t lf_den;

  size_t size; // number of set keys
  size_t capacity;
  size_t mask;

  void *entries;
  size_t *packed;

public:
  // instrumentation
  size_t lookups;
  size_t probes;

  HashTable(uint8_t tuple_size, uint8_t lf_num = DEFAULT_LF_NUM,
            uint8_t lf_den = DEFAULT_LF_DEN, size_t capacity = INITIAL_CAPACITY)
      : tuple_size(tuple_size),
        entry_size(tuple_size * sizeof(Key) + sizeof(Value)), lf_num(lf_num),
        lf_den(lf_den), size(0), capacity(capacity), mask(capacity - 1),
        lookups(0), probes(0) {
    entries = (void *)::operator new(capacity * entry_size);
    packed = new size_t[capacity];

    // memset to EMPTY_KEY
    for (size_t i = 0; i < capacity; i++) {
      Key *key = keyAt(i);
      key[0] = EMPTY_KEY; // this is sufficient
    }
  }

  ~HashTable() {
    ::operator delete(entries);
    delete[] packed;
  }

  void insert(const Key *key, Value value) {
    if ((size + 1) * lf_den >= lf_num * capacity) {
      resize();
    }
    size_t h = hash(key);
    while (keyAt(h)[0] != EMPTY_KEY) {
      h = (h + 1) & mask;
    }
    memcpy((char *)entries + h * entry_size, key, tuple_size * sizeof(Key));
    (*valueAt(h)) = value;
    packed[size++] = h;
  }

  Value *lookup(const Key *key) {
    if constexpr (Instrument) {
      ++lookups;
      ++probes;
    }
    size_t h = hash(key);
    Key *k = keyAt(h);
    while (k[0] != EMPTY_KEY) {
      if (keyEquals(key, k)) {
        return valueAt(h);
      }
      h = (h + 1) & mask;
      k = keyAt(h);
      if constexpr (Instrument) {
        ++probes;
      }
    }
    return nullptr;
  }

  // TODO: iter()

  size_t getSize() { return size; }

  size_t getCapacity() { return capacity; }

  size_t getMemoryUsage() { return capacity * entry_size; }

protected:
  [[gnu::always_inline]] Key *keyAt(size_t idx) {
    return (Key *)((char *)entries + idx * entry_size);
  }

  [[gnu::always_inline]] Key *keyAt(size_t idx, void *array) {
    return (Key *)((char *)array + idx * entry_size);
  }

  [[gnu::always_inline]] Value *valueAt(size_t idx) {
    return (Value *)((char *)entries + idx * entry_size +
                     tuple_size * sizeof(Key));
  }

  [[gnu::always_inline]] bool keyEquals(const Key *x, const Key *y) {
    for (size_t i = 0; i < tuple_size; i++) {
      if (x[i] != y[i]) {
        return false;
      }
    }
    return true;
  }

  void resize() {
    capacity <<= 1;
    mask = capacity - 1;

    void *new_entries = (void *)::operator new(capacity * entry_size);
    size_t *new_packed = new size_t[capacity];

    // memset to EMPTY_KEY
    for (size_t i = 0; i < capacity; i++) {
      Key *key = keyAt(i, new_entries);
      key[0] = EMPTY_KEY; // this is sufficient
    }

    // rehash
    for (size_t i = 0; i < size; ++i) {
      size_t idx = packed[i];
      Key *key = keyAt(idx);
      size_t h_new = hash(key);
      while (keyAt(h_new, new_entries)[0] != EMPTY_KEY) {
        h_new = (h_new + 1) & mask;
      }
      memcpy((char *)new_entries + h_new * entry_size,
             (char *)entries + idx * entry_size, entry_size);
      new_packed[i] = h_new;
    }

    ::operator delete(entries);
    delete[] packed;
    entries = new_entries;
    packed = new_packed;
  }

  size_t hash(const Key *key) {
    size_t h = 0;
    // TODO: idk man!
    for (size_t i = 0; i < tuple_size; i++) {
      h = h * 107 + key[i];
    }
    return h & mask;
  }
};

} // namespace space
