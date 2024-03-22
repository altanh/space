#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>

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
  static constexpr uint8_t DEFAULT_LF_DEN = 4;

  uint8_t tuple_size;
  const size_t entry_size;
  uint8_t lf_num;
  uint8_t lf_den;

  size_t size; // number of set keys
  size_t capacity;
  size_t mask;

  std::unique_ptr<Key[]> keys;
  std::unique_ptr<Value[]> values;
  std::unique_ptr<size_t[]> packed;

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
    keys = std::make_unique<Key[]>(capacity * tuple_size);
    values = std::make_unique<Value[]>(capacity);
    packed = std::make_unique<size_t[]>(capacity);

    for (size_t i = 0; i < capacity; i++) {
      keyAt(i)[0] = EMPTY_KEY; // this is sufficient
    }
  }

  // NB: we assume that the key is not already in the table
  void insert(const Key *key, Value value) {
    if ((size + 1) * lf_den >= lf_num * capacity) {
      resize();
    }
    size_t h = hash(key);
    while (keyAt(h)[0] != EMPTY_KEY) {
      h = (h + 1) & mask;
    }
    std::copy(key, key + tuple_size, keyAt(h));
    (*valueAt(h)) = value;
    packed[size++] = h;
  }

  void emplace(const Key *key, Value &&value) {
    if ((size + 1) * lf_den >= lf_num * capacity) {
      resize();
    }
    size_t h = hash(key);
    while (keyAt(h)[0] != EMPTY_KEY) {
      h = (h + 1) & mask;
    }
    std::copy(key, key + tuple_size, keyAt(h));
    (*valueAt(h)) = std::move(value);
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

  template <class KeyOut> void iter(size_t idx, KeyOut *out) {
    Key *k = keyAt(packed[idx]);
    // there might be a better way but this will do
    // if constexpr (std::is_same_v<Key, KeyOut>) {
    std::copy(k, k + tuple_size, out);
    // } else {
    //   for (size_t i = 0; i < tuple_size; i++) {
    //     out[i] = static_cast<KeyOut>(k[i]);
    //   }
    // }
  }

  template <class KeyOut, class ValueOut>
  void iter(size_t idx, KeyOut *out, ValueOut *val) {
    Key *k = keyAt(packed[idx]);
    // if constexpr (std::is_same_v<Key, KeyOut>) {
    std::copy(k, k + tuple_size, out);
    // } else {
    //   for (size_t i = 0; i < tuple_size; i++) {
    //     out[i] = static_cast<KeyOut>(k[i]);
    //   }
    // }
    *val = static_cast<ValueOut>(*valueAt(packed[idx]));
  }

  size_t getSize() { return size; }

  size_t getCapacity() { return capacity; }

  size_t getMemoryUsage() { return capacity * entry_size; }

protected:
  [[gnu::always_inline]] Key *keyAt(size_t idx) {
    return &keys[idx * tuple_size];
  }

  [[gnu::always_inline]] Key *keyAt(size_t idx, std::unique_ptr<Key[]> &array) {
    return &array[idx * tuple_size];
  }

  [[gnu::always_inline]] Value *valueAt(size_t idx) { return &values[idx]; }

  [[gnu::always_inline]] Value *valueAt(size_t idx,
                                        std::unique_ptr<Value[]> &array) {
    return &array[idx];
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

    std::unique_ptr<Key[]> new_keys(new Key[capacity * tuple_size]);
    std::unique_ptr<Value[]> new_values(new Value[capacity]);
    std::unique_ptr<size_t[]> new_packed(new size_t[capacity]);

    // memset to EMPTY_KEY
    for (size_t i = 0; i < capacity; i++) {
      Key *key = keyAt(i, new_keys);
      key[0] = EMPTY_KEY; // this is sufficient
    }

    // rehash
    for (size_t i = 0; i < size; ++i) {
      size_t idx = packed[i];
      Key *key = keyAt(idx);
      size_t h_new = hash(key);
      while (keyAt(h_new, new_keys)[0] != EMPTY_KEY) {
        h_new = (h_new + 1) & mask;
      }
      std::copy(key, key + tuple_size, keyAt(h_new, new_keys));
      (*valueAt(h_new, new_values)) = std::move(*valueAt(idx)); // sus?
      new_packed[i] = h_new;
    }

    // move new pointers
    keys = std::move(new_keys);
    values = std::move(new_values);
    packed = std::move(new_packed);
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

template <class Key, class Value, uint8_t TupleSize> class StaticHashTable {
  static constexpr Key EMPTY_KEY = std::numeric_limits<Key>::max();
  static constexpr size_t INITIAL_CAPACITY = 16;
  static constexpr uint8_t DEFAULT_LF_NUM = 1;
  static constexpr uint8_t DEFAULT_LF_DEN = 4;

  static constexpr size_t entry_size = TupleSize * sizeof(Key) + sizeof(Value);
  uint8_t lf_num;
  uint8_t lf_den;

  size_t size; // number of set keys
  size_t capacity;
  size_t mask;

  using Tuple = std::array<Key, TupleSize>;
  template <class T> using unique_ptr = std::unique_ptr<T>;

  unique_ptr<Tuple[]> keys;
  unique_ptr<Value[]> values;
  unique_ptr<size_t[]> packed;

public:
  StaticHashTable(uint8_t lf_num = DEFAULT_LF_NUM,
                  uint8_t lf_den = DEFAULT_LF_DEN,
                  size_t capacity = INITIAL_CAPACITY)
      : lf_num(lf_num), lf_den(lf_den), size(0), capacity(capacity),
        mask(capacity - 1) {
    keys = std::make_unique<Tuple[]>(capacity);
    values = std::make_unique<Value[]>(capacity);
    packed = std::make_unique<size_t[]>(capacity);

    for (size_t i = 0; i < capacity; i++) {
      keys[i][0] = EMPTY_KEY; // this is sufficient
    }
  }

  // NB: we assume that the key is not already in the table
  void insert(const Tuple &key, Value value) {
    if ((size + 1) * lf_den >= lf_num * capacity) {
      resize();
    }
    size_t h = hash(key);
    while (keys[h][0] != EMPTY_KEY) {
      h = (h + 1) & mask;
    }
    keys[h] = key;
    values[h] = value;
    packed[size++] = h;
  }

  void emplace(const Tuple &key, Value &&value) {
    if ((size + 1) * lf_den >= lf_num * capacity) {
      resize();
    }
    size_t h = hash(key);
    while (keys[h][0] != EMPTY_KEY) {
      h = (h + 1) & mask;
    }
    keys[h] = key;
    values[h] = std::move(value);
    packed[size++] = h;
  }

  Value *lookup(const Tuple &key) {
    size_t h = hash(key);
    while (keys[h][0] != EMPTY_KEY) {
      if (keyEquals(key, keys[h])) {
        return &values[h];
      }
      h = (h + 1) & mask;
    }
    return nullptr;
  }

  template <class KeyOut>
  void iter(size_t idx, std::array<KeyOut, TupleSize> &out) {
    if constexpr (std::is_same_v<Key, KeyOut>) {
      out = keys[packed[idx]];
    } else {
      std::copy(keys[packed[idx]].begin(), keys[packed[idx]].end(),
                out.begin());
    }
  }

  template <class KeyOut, class ValueOut>
  void iter(size_t idx, std::array<KeyOut, TupleSize> &out, ValueOut &val) {
    if constexpr (std::is_same_v<Key, KeyOut>) {
      out = keys[packed[idx]];
    } else {
      std::copy(keys[packed[idx]].begin(), keys[packed[idx]].end(),
                out.begin());
    }
    val = static_cast<ValueOut>(values[packed[idx]]);
  }

  Value *iter(size_t idx) { return &values[packed[idx]]; }

  // TODO: standardize naming
  size_t getSize() { return size; }

  size_t getCapacity() { return capacity; }

  size_t getMemoryUsage() { return capacity * entry_size; }

protected:
  bool keyEquals(const Tuple &x, const Tuple &y) {
    for (size_t i = 0; i < TupleSize; i++) {
      if (x[i] != y[i]) {
        return false;
      }
    }
    return true;
  }

  void resize() {
    capacity <<= 1;
    mask = capacity - 1;

    unique_ptr<Tuple[]> new_keys(new Tuple[capacity]);
    unique_ptr<Value[]> new_values(new Value[capacity]);
    unique_ptr<size_t[]> new_packed(new size_t[capacity]);

    // memset to EMPTY_KEY
    for (size_t i = 0; i < capacity; i++) {
      new_keys[i][0] = EMPTY_KEY; // this is sufficient
    }

    // rehash
    for (size_t i = 0; i < size; ++i) {
      size_t idx = packed[i];
      Tuple &key = keys[idx];
      size_t h_new = hash(key);
      while (new_keys[h_new][0] != EMPTY_KEY) {
        h_new = (h_new + 1) & mask;
      }
      new_keys[h_new] = key;
      new_values[h_new] = std::move(values[idx]);
      new_packed[i] = h_new;
    }

    // move new pointers
    keys = std::move(new_keys);
    values = std::move(new_values);
    packed = std::move(new_packed);
  }

  size_t hash(const Tuple &key) {
    size_t h = 0;
    // TODO: idk man!
    for (size_t i = 0; i < TupleSize; i++) {
      h = h * 107 + key[i];
    }
    return h & mask;
  }
};

} // namespace space
