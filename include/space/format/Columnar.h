#pragma once

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>
#include <vector>

namespace space {

namespace detail {

template <class NT> struct Values;
template <> struct Values<void> {};
template <class NT> struct Values {
  std::unique_ptr<NT[]> data;
  Values() : data(nullptr) {}
};

} // namespace detail

template <class IT, class NT> class Columnar {

  size_t rank;
  size_t size;
  // IT *dims;
  std::unique_ptr<IT[]> dims;
  detail::Values<NT> values;

public:
  static constexpr bool has_values = !std::is_void<NT>::value;

  Columnar(size_t size, size_t rank) {
    this->size = size;
    this->rank = rank;
    // dims = new IT[size * rank];
    dims = std::make_unique<IT[]>(size * rank);
    if constexpr (has_values) {
      values.data = std::make_unique<NT[]>(size);
    }
  }

  template <class NTIn, size_t Rank>
  static Columnar<IT, NT>
  fromCOO(const std::vector<std::pair<std::array<IT, Rank>, NTIn>> &coo) {
    size_t size = coo.size();
    size_t rank = Rank;
    Columnar<IT, NT> result(size, rank);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < rank; j++) {
        result.dim(j)[i] = coo[i].first[j];
      }
      if constexpr (has_values) {
        result.values.data[i] = coo[i].second;
      }
    }
    return result;
  }

  template <size_t Rank, class ValueInitializer>
  static Columnar<IT, NT> fromCOO(const std::vector<std::array<IT, Rank>> &coo,
                                  ValueInitializer init) {
    size_t size = coo.size();
    size_t rank = Rank;
    Columnar<IT, NT> result(size, rank);
    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < rank; j++) {
        result.dim(j)[i] = coo[i][j];
      }
      if constexpr (!std::is_void<NT>::value) {
        result.values.data[i] = init(coo[i]);
      }
    }
    return result;
  }

  IT *dim(size_t i) { return &dims[size * i]; }
  NT *val(size_t i) { return &values.data[i]; }

  size_t getSize() const { return size; }
  size_t getRank() const { return rank; }
};

} // namespace space
