#include <chrono>
#include <space/format/Columnar.h>
#include <space/format/Trie.h>

#include "Util.h"

// #define OLD

#ifdef OLD
#include <space/algorithm/FreeJoin_old.h>
#else
#include <space/algorithm/FreeJoin.h>
#endif

using namespace space;
using IT = uint64_t;
using NT = int64_t;

#ifdef OLD
void testJIT() {
  using std::array;
  using std::cout;
  using std::endl;
  using std::index_sequence;
  using std::tuple;
  using std::vector;

  // constexpr size_t N = 4000;
  constexpr double p = 0.1;

  vector<size_t> Ns = {3100};

  cout << "engine,n,m,triangles,time_s" << endl;

  for (auto N : Ns) {
    auto data = Columnar<IT, NT>(erdosRenyi<IT, NT>(N, p));
    const size_t M = data.getSize();

    cout << "fj," << N << "," << M << ",";

    // [R(x, y), S(y), T(x)], [S(z), T(z)]
    // using SR = index_sequence<2>;
    // using SS = index_sequence<1, 1>;
    // using ST = SS;
    // using PR = index_sequence<0, 1>;
    // using PS = PR;
    // using PT = PR;
    // using XY = index_sequence<0, 1>;
    // using YZ = index_sequence<1, 2>;
    // using XZ = index_sequence<0, 2>;
    // using N0 = index_sequence<0, 1, 2>;
    // using N1 = index_sequence<1, 2>;
    // using Nodes = jit::NodeList<N0, N1>;

    // [R(x), T(x)], [R(y), S(y)], [S(z), T(z)]
    // Schema
    using SR = index_sequence<1, 1>;
    using SS = SR;
    using ST = SR;
    // Data Permutation
    using PR = index_sequence<0, 1>;
    using PS = PR;
    using PT = PR;
    // Output Permutation
    using XY = index_sequence<0, 1>;
    using YZ = index_sequence<1, 2>;
    using XZ = index_sequence<0, 2>;
    // Nodes
    using N0 = index_sequence<0, 2>;
    using N1 = index_sequence<0, 1>;
    using N2 = index_sequence<1, 2>;
    using Nodes = jit::NodeList<N0, N1, N2>;

    StaticLazyTrie<IT, NT, true, true, SR, PR, XY> R(&data);
    StaticLazyTrie<IT, NT, true, true, SS, PS, YZ> S(&data);
    StaticLazyTrie<IT, NT, true, true, ST, PT, XZ> T(&data);
    tuple<decltype(R) *, decltype(S) *, decltype(T) *> rels(&R, &S, &T);

    size_t triangles = 0;

    // COUNT(*)
    auto count_tri = [&triangles](std::array<IT, 3> &idx, Value<NT> value) {
      triangles++;
    };

    // "CSR"-like output representation for masked-SpGEMM
    // vector<StaticHashTable<IT, NT, 1>> output(N);
    // auto mspgemm = [&output](array<IT, 3> &idx, Value<NT> value) {
    //   // project out y
    //   auto &row = output[idx[0]];
    //   auto col_val = row.lookup({idx[2]});
    //   if (col_val) {
    //     *col_val += value.some;
    //   } else {
    //     row.insert({idx[2]}, value.some);
    //   }
    // };

    // time it
    auto start = std::chrono::high_resolution_clock::now();
    array<IT, 3> tup;
    Value<NT> val;
    jit::FreeJoin<IT, NT, 3, decltype(rels), Nodes>::run(rels, tup, val,
                                                         count_tri);
    // reduce the matrix
    // for (size_t i = 0; i < N; i++) {
    //   auto &row = output[i];
    //   const size_t row_size = row.getSize();
    //   for (size_t j = 0; j < row_size; j++) {
    //     std::array<IT, 1> unused;
    //     NT v;
    //     row.iter(j, unused, v);
    //     triangles += v;
    //   }
    // }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << triangles << "," << elapsed_seconds.count() << std::endl;

    // GrB_Matrix A = util::toGraphBLAS<IT, NT>(data, false);

    // std::cout << "gb," << N << "," << M << ",";
    // GrB_Matrix C;
    // GrB_Matrix_new(&C, util::GraphBLASType<NT>::type(), N, N);
    // start = std::chrono::high_resolution_clock::now();
    // GrB_mxm(C, A, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_INT64, A, A, GrB_NULL);
    // int64_t gb_triangles;
    // GrB_Matrix_reduce_INT64(&gb_triangles, GrB_NULL, GrB_PLUS_MONOID_INT64,
    // C,
    //                         GrB_NULL);
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_seconds = end - start;
    // std::cout << gb_triangles << "," << elapsed_seconds.count() << std::endl;

    // // free A, C
    // GrB_Matrix_free(&A);
    // GrB_Matrix_free(&C);
  }
}
#else
void testJIT() {
  auto R_data = erdosRenyi2<IT, NT>(3100, 0.1);
  auto S_data = erdosRenyi2<IT, NT>(3100, 0.08, 1);
  auto T_data = erdosRenyi2<IT, NT>(3100, 0.3, 2);
  // auto data_T = data.permute({1, 0});

  // [R(x), T(x)], [R(y), S(y)], [S(z), T(z)]
  // auto R = jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<0, 1>, 2>(&data);
  // auto S = jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<1, 2>, 2>(&data_T);
  // auto T = jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<0, 2>, 2>(&data);

  std::array<IT, 3> tup;
  Value<NT> val;

  constexpr size_t nwarm = 2;
  constexpr size_t nrun = 10;

  std::vector<double> times;

  for (size_t i = 0; i < nwarm + nrun; ++i) {
    constexpr size_t R = 0;
    constexpr size_t S = 1;
    constexpr size_t T = 2;
    constexpr size_t X = 0;
    constexpr size_t Y = 1;
    constexpr size_t Z = 2;

    auto R_trie =
        jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<X, Y>, 2>(&R_data);
    auto S_trie =
        jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<Y, Z>, 2>(&S_data);
    auto T_trie =
        jit::LazyTrie<IT, NT, fn::List<1, 1>, fn::List<X, Z>, 2>(&T_data);

    // R.materialize();
    // S.materialize();
    // T.materialize();

    auto rels = std::make_tuple(&R_trie, &S_trie, &T_trie);

    constexpr size_t batch_size = 128;

    NT count = 0;
    auto start = std::chrono::high_resolution_clock::now();
    jit::FreeJoin<IT, NT, fn::List<X, Y, Z>,
                  jit::Nodes<fn::List<R, T>, fn::List<R, S>, fn::List<S, T>>,
                  batch_size>::run(rels, tup, val,
                                   [&count](auto t, auto v) { ++count; });
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "FreeJoin: " << count << " in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                     .count()
              << "ms" << std::endl;
    // if (count % 6 != 0) {
    //   std::cerr << "triangles should be divisible by 6" << std::endl;
    //   exit(1);
    // }
    if (i >= nwarm) {
      times.push_back(
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count());
    }
  }

  // compute mean and standard deviation
  double mean = 0;
  for (auto t : times) {
    mean += t;
  }
  mean /= nrun;

  double stddev = 0;
  for (auto t : times) {
    stddev += (t - mean) * (t - mean);
  }
  stddev = std::sqrt(stddev / nrun);

  std::cout << "mean: " << mean << "ms" << std::endl;
  std::cout << "stddev: " << stddev << "ms" << std::endl;
}
#endif

int main() {
  testJIT();

  return 0;
}