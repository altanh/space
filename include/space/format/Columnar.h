#pragma once

#include <type_traits>

template <class IT, class NT> class Columnar {
  IT nnz;
  IT ndim;
  IT *dims;
  NT *values;

public:
  static constexpr bool has_values = !std::is_void<NT>::value;

  Columnar(IT nnz, IT ndim) {
    this->nnz = nnz;
    this->ndim = ndim;
    dims = new IT[nnz * ndim];
    if constexpr (std::is_void<NT>::value) {
      values = nullptr;
    } else {
      values = new NT[nnz];
    }
  }

  IT *dim(IT i) { return dims + nnz * i; }
};
