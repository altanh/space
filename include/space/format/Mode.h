#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <limits>
#include <optional>
#include <tuple>
#include <vector>

#include "HashTable.h"

namespace space {

constexpr size_t Unbounded = std::numeric_limits<size_t>::max();

template <class D, class V, size_t U, class Impl> class ModeIterator {
public:
  virtual bool hasNext() = 0;
  virtual Impl &operator++() = 0;
  virtual std::pair<const D &, V &> operator*() = 0;
  virtual const D &key() = 0;
  virtual V &value() = 0;
};

template <class D, class V, size_t U, class _Iterator> class Mode {
public:
  using Iterator = _Iterator;

  virtual V *insert(const D &key, V value) = 0;
  virtual V *emplace(const D &key, V &&value) = 0;
  virtual V *lookup(const D &key) = 0;
  virtual void finalize() = 0;

  virtual Iterator iterator() = 0;
};

// Helper to compute the product of an index sequence.
template <size_t... Is> struct Product;
template <size_t I, size_t... Is> struct Product<I, Is...> {
  static constexpr size_t value = I * Product<Is...>::value;
};

template <> struct Product<> {
  static constexpr size_t value = 1;
};

// Fused mode.
template <class D, class V, class _Iterator,
          template <typename, typename, size_t, typename> class M, size_t... Us>
class FusedMode : public Mode<std::array<D, sizeof...(Us)>, V,
                              Product<Us...>::value, _Iterator> {
  constexpr static size_t U = Product<Us...>::value;
  constexpr static size_t K = sizeof...(Us);

  M<D, V, U, _Iterator> mode_;

  // TODO: make this work somehow
};

template <class D, class V, size_t U> class DenseMode;

template <class D, class V, size_t U>
class DenseModeIterator
    : public ModeIterator<D, V, U, DenseModeIterator<D, V, U>> {
  DenseMode<D, V, U> &mode_;
  size_t idx_;

public:
  DenseModeIterator(DenseMode<D, V, U> &mode) : mode_(mode), idx_(0) {}

  bool hasNext() override { return idx_ < mode_.nonzeros_.size(); }

  DenseModeIterator &operator++() override {
    idx_++;
    return *this;
  }

  std::pair<const D &, V &> operator*() override {
    return {mode_.nonzeros_[idx_], *mode_.values_[mode_.nonzeros_[idx_]]};
  }

  const D &key() override { return mode_.nonzeros_[idx_]; }

  V &value() override { return *mode_.values_[mode_.nonzeros_[idx_]]; }
};

template <class D, class V, size_t U>
class DenseMode : public Mode<D, V, U, DenseModeIterator<D, V, U>> {
  friend class DenseModeIterator<D, V, U>;

  static_assert(U != Unbounded, "DenseMode must have a bounded capacity");
  // Not sure how to avoid optional here without sentinel value for V.
  // Alternative is to use a bitset to track non-empty values.
  std::vector<std::optional<V>> values_;
  std::vector<D> nonzeros_;

public:
  using Iterator = DenseModeIterator<D, V, U>;

  DenseMode() : values_(U, std::nullopt), nonzeros_() {}

  V *insert(const D &key, V value) override {
    values_[key] = value;
    nonzeros_.push_back(key);
    return values_[key].operator->();
  }

  V *emplace(const D &key, V &&value) override {
    values_[key] = std::move(value);
    nonzeros_.push_back(key);
    return values_[key].operator->();
  }

  V *lookup(const D &key) override {
    if (values_[key]) {
      return values_[key].operator->();
    }
    return nullptr;
  }

  void finalize() override {}

  Iterator iterator() override { return Iterator(*this); }
};

// Dense mode alternative using a bitset.
template <class D, class V, size_t U> class BitsetDenseMode;

template <class D, class V, size_t U>
class BitsetDenseModeIterator
    : public ModeIterator<D, V, U, BitsetDenseModeIterator<D, V, U>> {
  BitsetDenseMode<D, V, U> &mode_;
  size_t idx_;

public:
  BitsetDenseModeIterator(BitsetDenseMode<D, V, U> &mode)
      : mode_(mode), idx_(0) {}

  bool hasNext() override { return idx_ < mode_.nonzeros_.size(); }

  BitsetDenseModeIterator &operator++() override {
    idx_++;
    return *this;
  }

  std::pair<const D &, V &> operator*() override {
    return {mode_.nonzeros_[idx_], mode_.values_[mode_.nonzeros_[idx_]]};
  }

  const D &key() override { return mode_.nonzeros_[idx_]; }

  V &value() override { return mode_.values_[mode_.nonzeros_[idx_]]; }
};

template <class D, class V, size_t U>
class BitsetDenseMode : public Mode<D, V, U, BitsetDenseModeIterator<D, V, U>> {
  friend class BitsetDenseModeIterator<D, V, U>;

  static_assert(U != Unbounded, "BitsetDenseMode must have a bounded capacity");
  std::vector<V> values_;
  std::bitset<U> bitset_;
  std::vector<D> nonzeros_;

public:
  using Iterator = BitsetDenseModeIterator<D, V, U>;

  BitsetDenseMode() : values_(U), bitset_(), nonzeros_() {}

  V *insert(const D &key, V value) override {
    values_[key] = value;
    bitset_[key] = true;
    nonzeros_.push_back(key);
    return &values_[key];
  }

  V *emplace(const D &key, V &&value) override {
    values_[key] = std::move(value);
    bitset_[key] = true;
    nonzeros_.push_back(key);
    return &values_[key];
  }

  V *lookup(const D &key) override {
    if (bitset_[key]) {
      return &values_[key];
    }
    return nullptr;
  }

  void finalize() override {}

  Iterator iterator() override { return Iterator(*this); }
};

// Sorted mode.
template <class D, class V, size_t U> class SortedMode;

template <class D, class V, size_t U>
class SortedModeIterator
    : public ModeIterator<D, V, U, SortedModeIterator<D, V, U>> {
  SortedMode<D, V, U> &mode_;
  size_t idx_;

public:
  SortedModeIterator(SortedMode<D, V, U> &mode) : mode_(mode), idx_(0) {}

  bool hasNext() override { return idx_ < mode_.entries_.size(); }

  SortedModeIterator &operator++() override {
    idx_++;
    return *this;
  }

  std::pair<const D &, V &> operator*() override {
    auto &entry = mode_.entries_[idx_];
    return {std::get<0>(entry), std::get<1>(entry)};
  }

  const D &key() override { return std::get<0>(mode_.entries_[idx_]); }

  V &value() override { return std::get<1>(mode_.entries_[idx_]); }
};

template <class D, class V, size_t U>
class SortedMode : public Mode<D, V, U, SortedModeIterator<D, V, U>> {
public:
  using Entry = std::tuple<D, V>;

private:
  friend class SortedModeIterator<D, V, U>;

  bool finalized_;
  size_t idx_;
  std::vector<Entry> entries_;

public:
  using Iterator = SortedModeIterator<D, V, U>;

  SortedMode() : finalized_(false), idx_(0), entries_() {}

  V *insert(const D &key, V value) override {
    assert(!finalized_);
    entries_.emplace_back(key, value);
    return &std::get<1>(entries_.back());
  }

  V *emplace(const D &key, V &&value) override {
    assert(!finalized_);
    entries_.emplace_back(key, std::move(value));
    return &std::get<1>(entries_.back());
  }

  V *lookup(const D &key) override {
    assert(finalized_);
    return lookupGallop(key);
  }

  void finalize() override {
    if (!finalized_) {
      std::sort(entries_.begin(), entries_.end());
      finalized_ = true;
    }
  }

  Iterator iterator() override { return Iterator(*this); }

protected:
  Entry *prepareGallop(const D &key, bool *forward) {
    if (idx_ == entries_.size()) {
      *forward = false;
    } else {
      auto &cur_key = std::get<0>(entries_[idx_]);
      if (cur_key == key) {
        return &entries_[idx_++];
      } else if (cur_key < key) {
        *forward = true;
      } else {
        *forward = false;
      }
    }
    return nullptr;
  }

  V *lookupGallop(const D &key) {
    if (entries_.empty()) {
      return nullptr;
    }
    bool forward;
    auto entry = prepareGallop(key, &forward);
    if (entry) {
      return &std::get<1>(*entry);
    }
    if (forward) {
      // forward galloping search
      // NB: prepareGallop guarantees idx_ < entries_.size(), so this is safe
      const size_t offset_max = entries_.size() - idx_ - 1;
      size_t offset = 1;
      while (offset < offset_max &&
             std::get<0>(entries_[idx_ + offset]) < key) {
        offset *= 2;
      };
      // increment idx_ since we already checked idx_ in prepareGallop
      idx_++;
      // binary search from idx_ + offset / 2 to idx_ + offset
      size_t lo = idx_ + offset / 2;
      size_t hi = idx_ + std::min(offset, offset_max);
      while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        idx_ = mid + 1;
        entry = &entries_[mid];
        if (std::get<0>(*entry) == key) {
          return &std::get<1>(*entry);
        } else if (std::get<0>(*entry) < key) {
          lo = mid + 1;
        } else {
          hi = mid;
        }
      }
    } else {
      // backward galloping search
      size_t offset = 1;
      while (offset <= idx_ && std::get<0>(entries_[idx_ - offset]) > key) {
        offset *= 2;
      }
      // binary search from idx_ - offset to idx_ - offset / 2
      size_t lo = idx_ - std::min(offset, idx_);
      size_t hi = idx_ - offset / 2;
      while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        idx_ = mid + 1;
        entry = &entries_[mid];
        if (std::get<0>(*entry) == key) {
          return &std::get<1>(*entry);
        } else if (std::get<0>(*entry) < key) {
          lo = mid + 1;
        } else {
          hi = mid;
        }
      }
    }
    return nullptr;
  }

  Entry *prepareBinary(const D &key, size_t *lo, size_t *hi) {
    if (idx_ == entries_.size()) {
      *lo = 0;
      *hi = entries_.size();
    } else {
      auto &cur_key = std::get<0>(entries_[idx_]);
      if (cur_key == key) {
        return &entries_[idx_++];
      } else if (cur_key < key) {
        *lo = idx_ + 1;
        *hi = entries_.size();
      } else {
        *lo = 0;
        *hi = idx_;
      }
    }
    return nullptr;
  }

  V *lookupBinary(const D &key) {
    if (entries_.empty()) {
      return nullptr;
    }
    size_t lo, hi;
    auto entry = prepareBinary(key, &lo, &hi);
    if (entry) {
      return &std::get<1>(*entry);
    }
    while (lo < hi) {
      size_t mid = lo + (hi - lo) / 2;
      idx_ = mid + 1;
      entry = &entries_[mid];
      if (std::get<0>(*entry) == key) {
        return &std::get<1>(*entry);
      } else if (std::get<0>(*entry) < key) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    return nullptr;
  }
};

// Hash mode.
template <class D, class V, size_t U> class HashMode;
template <class D, class V, size_t U> class HashModeIterator;

template <class D, class V, size_t U, size_t K>
class HashModeIterator<std::array<D, K>, V, U>
    : public ModeIterator<std::array<D, K>, V, U,
                          HashModeIterator<std::array<D, K>, V, U>> {
  HashMode<std::array<D, K>, V, U> &mode_;
  size_t idx_;

public:
  HashModeIterator(HashMode<std::array<D, K>, V, U> &mode)
      : mode_(mode), idx_(0) {}

  bool hasNext() override { return idx_ < mode_.table_.getSize(); }

  HashModeIterator &operator++() override {
    idx_++;
    return *this;
  }

  std::pair<const std::array<D, K> &, V &> operator*() override {
    auto entry = mode_.table_.iterEntry(idx_);
    return {entry.first, entry.second};
  }

  const std::array<D, K> &key() override {
    return mode_.table_.iterEntry(idx_).first;
  }

  V &value() override { return mode_.table_.iterEntry(idx_).second; }
};

template <class D, class V, size_t U>
class HashModeIterator
    : public ModeIterator<D, V, U, HashModeIterator<D, V, U>> {
  HashMode<D, V, U> &mode_;
  size_t idx_;

public:
  HashModeIterator(HashMode<D, V, U> &mode) : mode_(mode), idx_(0) {}

  bool hasNext() override { return idx_ < mode_.table_.getSize(); }

  HashModeIterator &operator++() override {
    idx_++;
    return *this;
  }

  std::pair<const D &, V &> operator*() override {
    auto entry = mode_.table_.iterEntry(idx_);
    // unwrap the length 1 array key
    return {std::get<0>(entry.first), entry.second};
  }

  const D &key() override {
    return std::get<0>(mode_.table_.iterEntry(idx_).first);
  }

  V &value() override { return mode_.table_.iterEntry(idx_).second; }
};

// Tuple specialization
template <class D, class V, size_t U, size_t K>
class HashMode<std::array<D, K>, V, U>
    : public Mode<std::array<D, K>, V, U,
                  HashModeIterator<std::array<D, K>, V, U>> {
  friend class HashModeIterator<std::array<D, K>, V, U>;

  StaticHashTable<D, V, K> table_;

public:
  using Iterator = HashModeIterator<std::array<D, K>, V, U>;
  using Key = std::array<D, K>;

  HashMode() : table_() {}

  V *insert(const Key &key, V value) override {
    return table_.insert(key, value);
  }

  V *emplace(const Key &key, V &&value) override {
    return table_.emplace(key, std::move(value));
  }

  V *lookup(const Key &key) override { return table_.lookup(key); }

  void finalize() override {}

  Iterator iterator() override { return Iterator(*this); }
};

template <class D, class V, size_t U>
class HashMode : public Mode<D, V, U, HashModeIterator<D, V, U>> {
  friend class HashModeIterator<D, V, U>;

  StaticHashTable<D, V, 1> table_;

public:
  using Iterator = HashModeIterator<D, V, U>;
  using Key = D;

  HashMode() : table_() {}

  V *insert(const Key &key, V value) override {
    return table_.insert({key}, value);
  }

  V *emplace(const Key &key, V &&value) override {
    return table_.emplace({key}, std::move(value));
  }

  V *lookup(const Key &key) override { return table_.lookup({key}); }

  void finalize() override {}

  Iterator iterator() override { return Iterator(*this); }
};

} // namespace space
