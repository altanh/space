#pragma once

#include <GraphBLAS.h>
#include <memory>
#include <numeric>
#include <space/format/Columnar.h>

namespace space {

namespace util {

template <typename NT> struct GraphBLASType;

template <> struct GraphBLASType<float> {
  static GrB_Type type() { return GrB_FP32; }

  // specialize GrB_Matrix_import for COO
  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, float *X,
                         GrB_Index nvals) {
    return GrB_Matrix_import_FP32(A, type(), nrows, ncols, I_, J, X, nvals,
                                  nvals, nvals, GrB_COO_FORMAT);
  }
};

template <> struct GraphBLASType<double> {
  static GrB_Type type() { return GrB_FP64; }

  // specialize GrB_Matrix_import for COO
  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, double *X,
                         GrB_Index nvals) {
    return GrB_Matrix_import_FP64(A, type(), nrows, ncols, I_, J, X, nvals,
                                  nvals, nvals, GrB_COO_FORMAT);
  }
};

template <> struct GraphBLASType<int> {
  static GrB_Type type() { return GrB_INT32; }

  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, int *X, GrB_Index nvals) {
    return GrB_Matrix_import_INT32(A, type(), nrows, ncols, I_, J, X, nvals,
                                   nvals, nvals, GrB_COO_FORMAT);
  }
};

template <> struct GraphBLASType<unsigned int> {
  static GrB_Type type() { return GrB_UINT32; }

  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, unsigned int *X,
                         GrB_Index nvals) {
    return GrB_Matrix_import_UINT32(A, type(), nrows, ncols, I_, J, X, nvals,
                                    nvals, nvals, GrB_COO_FORMAT);
  }
};

template <> struct GraphBLASType<long> {
  static GrB_Type type() { return GrB_INT64; }

  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, long *X,
                         GrB_Index nvals) {
    return GrB_Matrix_import_INT64(A, type(), nrows, ncols, I_, J, X, nvals,
                                   nvals, nvals, GrB_COO_FORMAT);
  }
};

template <> struct GraphBLASType<unsigned long> {
  static GrB_Type type() { return GrB_UINT64; }

  static GrB_Info import(GrB_Matrix *A, GrB_Index nrows, GrB_Index ncols,
                         GrB_Index *I_, GrB_Index *J, unsigned long *X,
                         GrB_Index nvals) {
    return GrB_Matrix_import_UINT64(A, type(), nrows, ncols, I_, J, X, nvals,
                                    nvals, nvals, GrB_COO_FORMAT);
  }
};

template <typename IT, typename NT>
GrB_Matrix toGraphBLAS(const Columnar<IT, NT> &columnar,
                       bool transpose = false) {
  if (columnar.getRank() != 2) {
    throw std::invalid_argument("Columnar matrix must have rank 2");
  }

  const IT row_dim = transpose ? 1 : 0;
  const IT col_dim = transpose ? 0 : 1;
  const IT nnz = columnar.getSize();

  const IT *row_indices = columnar.dim(row_dim);
  const IT *col_indices = columnar.dim(col_dim);

  // assuming 0 based indexing
  const IT nrows = std::accumulate(row_indices, row_indices + nnz, 0,
                                   [](IT a, IT b) { return std::max(a, b); }) +
                   1;
  const IT ncols = std::accumulate(col_indices, col_indices + nnz, 0,
                                   [](IT a, IT b) { return std::max(a, b); }) +
                   1;

  std::unique_ptr<GrB_Index[]> r = std::make_unique<GrB_Index[]>(nnz);
  std::unique_ptr<GrB_Index[]> c = std::make_unique<GrB_Index[]>(nnz);
  std::copy(row_indices, row_indices + nnz, r.get());
  std::copy(col_indices, col_indices + nnz, c.get());

  GrB_Matrix A;
  GraphBLASType<NT>::import(&A, nrows, ncols, c.get(), r.get(), columnar.val(0),
                            nnz);

  return A;
}

} // namespace util

} // namespace space