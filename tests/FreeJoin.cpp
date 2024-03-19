#include <array>
#include <chrono>
#include <cmath>
#include <space/algorithm/FreeJoin.h>
#include <space/util/GB.h>
#include <vector>

#include <GraphBLAS.h>

#include "Util.h"

using namespace space;
using std::array;
using std::vector;
using IT = uint64_t;
using NT = double;

#include <regex>
#include <sstream>

FreeJoinPlan parsePlan(std::string plan) {
  // Grammar:
  // plan -> node (',' node)*
  // node -> '[' atom (',' atom)* ']'
  // atom -> var '(' var (',' var)* ')'
  // var -> [a-zA-Z_][a-zA-Z0-9_]*
  FreeJoinPlan result;

  std::regex node_pattern("\\[([^\\]]+)\\]");
  std::regex atom_pattern("([a-zA-Z_][a-zA-Z0-9_]*)\\(([^\\)]+)\\)");
  std::regex var_pattern("[a-zA-Z_][a-zA-Z0-9_]*");

  for (std::sregex_iterator it(plan.begin(), plan.end(), node_pattern);
       it != std::sregex_iterator(); ++it) {
    std::smatch node_match = *it;
    std::string node = node_match[1];

    std::vector<Atom> atoms;
    for (std::sregex_iterator jt(node.begin(), node.end(), atom_pattern);
         jt != std::sregex_iterator(); ++jt) {
      std::smatch atom_match = *jt;
      std::string atom_name = atom_match[1];
      std::string vars = atom_match[2];

      std::vector<std::string> var_names;
      for (std::sregex_iterator kt(vars.begin(), vars.end(), var_pattern);
           kt != std::sregex_iterator(); ++kt) {
        std::smatch var_match = *kt;
        var_names.push_back(var_match[0]);
      }

      atoms.push_back(Atom(atom_name, var_names));
    }

    result.push_back(atoms);
  }

  // print the parsed plan
  // std::cout << "Parsed plan:" << std::endl;
  // for (const auto &node : result) {
  //   std::cout << "[ ";
  //   for (const auto &atom : node) {
  //     std::cout << atom.relation << "(";
  //     for (size_t i = 0; i < atom.vars.size(); i++) {
  //       std::cout << atom.vars[i];
  //       if (i + 1 < atom.vars.size()) {
  //         std::cout << " ";
  //       }
  //     }
  //     std::cout << ") ";
  //   }
  //   std::cout << "] ";
  // }
  // std::cout << std::endl;

  return result;
}

void testER() {
  // constexpr size_t N = 4000;
  constexpr double p = 0.1;

  std::vector<size_t> Ns;
  for (size_t i = 100; i <= 3100; i += 500) {
    Ns.push_back(i);
  }

  auto plan = parsePlan("[R(x, y), S(y), T(x)], [S(z), T(z)]");
  // auto plan = parsePlan("[R(x), T(x)], [R(y), S(y)], [S(z), T(z)]");
  // auto plan = parsePlan("[R(x, y), S(y)], [S(z), T(z)]");

  std::cout << "engine,n,m,triangles,time_s" << std::endl;

  GrB_init(GrB_NONBLOCKING);
  GxB_Context_set(GxB_CONTEXT_WORLD, GxB_CONTEXT_NTHREADS, 1);

  for (auto N : Ns) {
    auto data = Columnar<IT, NT>(erdosRenyi<IT, NT>(N, p));
    const size_t M = data.getSize();

    // std::cout << "Edges: " << M << std::endl;
    // std::cout << "Triangle upper bound: " << (size_t)std::pow(M, 1.5)
    //           << std::endl;
    std::cout << "fj," << N << "," << M << ",";

    auto R = Relation<IT, NT>{"R", {"x", "y"}, &data};
    auto S = Relation<IT, NT>{"S", {"y", "z"}, &data};
    auto T = Relation<IT, NT>{"T", {"z", "x"}, &data};

    std::vector<Relation<IT, NT>> all_relations{R, S, T};
    size_t triangles = 0;

    // time it
    auto start = std::chrono::high_resolution_clock::now();
    FreeJoin(all_relations, plan,
             [&triangles](IT *idx, Value<NT> value) { triangles++; });
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << triangles << "," << elapsed_seconds.count() << std::endl;

    if (triangles % 6 != 0) {
      std::cerr << "Error: " << triangles << " is not divisible by 6"
                << std::endl;
      exit(1);
    }
    // std::cout << "Triangles: " << triangles / 6 << std::endl;

    GrB_Matrix A = util::toGraphBLAS<IT, NT>(data, false);

    std::cout << "gb," << N << "," << M << ",";
    GrB_Matrix C;
    GrB_Matrix_new(&C, util::GraphBLASType<NT>::type(), N, N);
    start = std::chrono::high_resolution_clock::now();
    GrB_mxm(C, A, GrB_NULL, GrB_PLUS_TIMES_SEMIRING_FP64, A, A, GrB_NULL);
    double gb_triangles;
    GrB_Matrix_reduce_FP64(&gb_triangles, GrB_NULL, GrB_PLUS_MONOID_FP64, C,
                           GrB_NULL);
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    std::cout << gb_triangles << "," << elapsed_seconds.count() << std::endl;

    // free A, C
    GrB_Matrix_free(&A);
    GrB_Matrix_free(&C);
  }
}

void testWorstCase() {
  std::vector<size_t> Ns;
  for (size_t i = 100; i <= 20000; i += 100) {
    Ns.push_back(i);
  }

  auto plan = parsePlan("[R(x, y), S(y), T(x)], [S(z), T(z)]");

  std::cout << "in,out,time_ms" << std::endl;

  for (auto N : Ns) {
    auto [Rd, Sd, Td] = triangleWorstCase<IT, NT>(N);
    auto R = Relation<IT, NT>{"R", {"x", "y"}, &Rd};
    auto S = Relation<IT, NT>{"S", {"y", "z"}, &Sd};
    auto T = Relation<IT, NT>{"T", {"z", "x"}, &Td};
    std::vector<Relation<IT, NT>> all_relations{R, S, T};
    size_t triangles = 0;

    // time it
    auto start = std::chrono::high_resolution_clock::now();
    FreeJoin(all_relations, plan,
             [&triangles](IT *idx, Value<NT> value) { triangles++; });
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << N << "," << triangles << "," << elapsed_seconds.count()
              << std::endl;

    // if (triangles % 6 != 0) {
    //   std::cerr << "Error: " << triangles << " is not divisible by 6"
    //             << std::endl;
    //   exit(1);
    // }
    // std::cout << "Triangles: " << triangles / 6 << std::endl;
  }
}

void testJIT() {
  {
    std::cout << "[R(x,y), S(y)], [S(z)]" << std::endl;
    constexpr int VO = 3;
    auto data = Columnar<IT, NT>(erdosRenyi<IT, NT>(4, 1.0));

    StaticLazyTrie<IT, NT, true, true, std::index_sequence<2>,
                   std::index_sequence<0, 1>, std::index_sequence<0, 1>>
        R(&data);

    static_assert(decltype(R)::OutPerm::size() == 2, "R perm size");

    StaticLazyTrie<IT, NT, true, true, std::index_sequence<1, 1>,
                   std::index_sequence<1, 0>, std::index_sequence<1, 2>>
        S(&data);

    std::tuple<decltype(R) *, decltype(S) *> rels(&R, &S);
    std::array<IT, 3> tup;
    Value<NT> val;

    auto output = [](std::array<IT, 3> &tup, Value<NT> &val) {
      std::cout << "output: " << tup[0] << ", " << tup[1] << ", " << tup[2]
                << " -> " << val.some << std::endl;
    };
    jit::FreeJoin<IT, NT, VO, decltype(rels),
                  jit::NodeList<std::index_sequence<0, 1>,
                                std::index_sequence<1>>>::run(rels, tup, val,
                                                              output);
  }

  // {
  //   std::cout << "[R(x,y)], [S(z)]" << std::endl;
  //   constexpr int VO = 3;
  //   constexpr int VP = 0;
  //   auto data = Columnar<IT, NT>(erdosRenyi<IT, NT>(3, 1.0));

  //   StaticLazyTrie<IT, NT, true, true, std::index_sequence<2>,
  //                  std::index_sequence<0, 1>>
  //       R(&data);

  //   StaticLazyTrie<IT, NT, true, true, std::index_sequence<1>,
  //                  std::index_sequence<0>>
  //       S(&data);

  //   std::tuple<decltype(R) *, decltype(S) *> rels(&R, &S);
  //   std::array<IT, 3> tup;
  //   Value<NT> val;

  //   jit::FreeJoin<IT, NT, VO, VP, decltype(rels),
  //                 jit::NodeList<std::index_sequence<0>,
  //                               std::index_sequence<1>>>::run(rels, tup,
  //                               val);
  // }

  // {
  //   std::cout << "[R(x,y), S(x,y)]" << std::endl;
  //   constexpr int VO = 2;
  //   auto data1 = Columnar<IT, NT>(erdosRenyi<IT, NT>(3, 1.0));
  //   auto data0 = Columnar<IT, NT>(erdosRenyi<IT, NT>(10, 1.0));

  //   StaticLazyTrie<IT, NT, true, true, std::index_sequence<2>,
  //                  std::index_sequence<0, 1>, std::index_sequence<0, 1>>
  //       R(&data0);

  //   StaticLazyTrie<IT, NT, true, true, std::index_sequence<2>,
  //                  std::index_sequence<0, 1>, std::index_sequence<0, 1>>
  //       S(&data1);

  //   std::tuple<decltype(R) *, decltype(S) *> rels(&R, &S);
  //   std::array<IT, 2> tup;
  //   Value<NT> val;

  //   jit::FreeJoin<IT, NT, VO, decltype(rels),
  //                 jit::NodeList<std::index_sequence<0, 1>>>::run(rels, tup,
  //                                                                val);
  // }
}

int main(int argc, char **argv) {
  // testER();
  // testWorstCase();
  testJIT();
  return 0;
}