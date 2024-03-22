#pragma once

#include <array>
#include <iostream>
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

namespace jit {

namespace detail {

template <class NT> struct Values;
template <> struct Values<void> {};
template <class NT> struct Values {
  std::shared_ptr<std::vector<NT>> data;
  Values() : data(nullptr) {}
};

} // namespace detail

template <class IT, class NT, size_t Rank> class Columnar {
  using Values = std::vector<NT>;
  using Column = std::vector<IT>;
  using ColumnPtr = std::shared_ptr<Column>;

  size_t size_;
  std::array<ColumnPtr, Rank> cols_;
  detail::Values<NT> vals_;

public:
  Columnar() : size_(0) {}

  // Copy constructor
  Columnar(const Columnar &other) { *this = other; }

  Columnar(Columnar &&other) {
    size_ = other.size_;
    other.size_ = 0;
    cols_ = std::move(other.cols_);
    if constexpr (hasValues()) {
      vals_.data = std::move(other.vals_.data);
    }
  }

  Columnar(size_t size) : size_(size) {
    for (size_t col = 0; col < Rank; ++col) {
      cols_[col] = std::make_shared<Column>(size);
    }
    if constexpr (hasValues()) {
      vals_.data = std::make_shared<Values>(size);
    }
  }

  // Copy assignment
  Columnar &operator=(const Columnar &other) {
    if (this != &other) {
      size_ = other.size_;
      cols_ = other.cols_;
      if constexpr (hasValues()) {
        vals_.data = other.vals_.data;
      }
    }
    return *this;
  }

  // Move assignment
  Columnar &operator=(Columnar &&other) {
    if (this != &other) {
      size_ = other.size_;
      other.size_ = 0;
      cols_ = std::move(other.cols_);
      if constexpr (hasValues()) {
        vals_.data = std::move(other.vals_.data);
      }
    }
    return *this;
  }

  template <class ValueInitializer>
  static Columnar fromCOO(const std::vector<std::array<IT, Rank>> &coo,
                          ValueInitializer &&init) {
    Columnar result(coo.size());
    for (size_t i = 0; i < coo.size(); ++i) {
      auto r = result.row(i);
      for (size_t col = 0; col < Rank; ++col) {
        r[col] = coo[i][col];
      }
      if constexpr (result.hasValues()) {
        r.value() = init(i, coo[i]);
      }
    }
    return result;
  }

  Column &column(size_t col) { return *cols_[col]; }

  Column &column(size_t col) const { return *cols_[col]; }

  NT &value(size_t index) { return (*vals_.data)[index]; }

  NT &value(size_t index) const { return (*vals_.data)[index]; }

  class RowView {
    size_t index;
    Columnar *self;

  public:
    RowView(size_t i, Columnar *s) : index(i), self(s) {}

    IT &operator[](size_t col) { return self->column(col)[index]; }
    IT &operator[](size_t col) const { return self->column(col)[index]; }
    NT &value() { return self->value(index); }
    NT &value() const { return self->value(index); }
  };

  RowView row(size_t index) { return {index, this}; }

  RowView row(size_t index) const { return {index, this}; }

  Columnar permute(std::array<size_t, Rank> perm) {
    Columnar result;
    result.size_ = size_;
    for (size_t col = 0; col < rank(); ++col) {
      result.cols_[perm[col]] = cols_[col];
    }
    result.vals_ = vals_;
    return result;
  }

  void print() const {
    for (size_t index = 0; index < size_; ++index) {
      std::cout << "(" << (*cols_[0])[index];
      for (size_t col = 1; col < rank(); ++col) {
        std::cout << "," << (*cols_[col])[index];
      }
      std::cout << ")";
      if constexpr (hasValues()) {
        std::cout << " -> " << (*vals_.data)[index];
      }
      std::cout << std::endl;
    }
  }

  static constexpr bool hasValues() { return !std::is_void_v<NT>; }

  static constexpr size_t rank() { return Rank; }

  size_t size() const { return size_; }
};

} // namespace jit

// TODO: consider dims as pointer to pointer
template <class IT, class NT> class Columnar {
  size_t rank;
  size_t size;
  std::unique_ptr<IT[]> dims;
  detail::Values<NT> values;

public:
  static constexpr bool has_values = !std::is_void<NT>::value;

  // Columnar(const Columnar<IT, NT> &other, std::unique_ptr<size_t[]>
  // permutation)
  //     : rank(other.rank), size(other.size),
  //     permutation(std::move(permutation)),
  //       dims(other.dims), values(other.values) {}

  Columnar(size_t size, size_t rank) {
    this->size = size;
    this->rank = rank;
    dims = std::make_unique<IT[]>(size * rank);
    // dims = std::shared_ptr<IT>(new IT[size * rank], [](IT *p) { delete[] p;
    // });
    if constexpr (has_values) {
      values.data = std::make_unique<NT[]>(size);
      // values.data =
      //     std::shared_ptr<NT>(new NT[size], [](NT *p) { delete[] p; });
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
  IT *dim(size_t i) const { return &dims[size * i]; }

  NT *val(size_t i) { return &values.data[i]; }
  NT *val(size_t i) const { return &values.data[i]; }

  size_t getSize() const { return size; }
  size_t getRank() const { return rank; }
};

} // namespace space
