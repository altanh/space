#pragma once

#include <array>
#include <random>
#include <space/format/Columnar.h>

template <class IT, class NT>
space::Columnar<IT, NT> erdosRenyi(size_t n, double p, size_t seed = 0) {
  std::mt19937 gen(seed);
  std::bernoulli_distribution d(p);

  std::vector<std::array<IT, 2>> coo;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (d(gen)) {
        coo.push_back({i, j});
        coo.push_back({j, i});
      }
    }
  }

  return space::Columnar<IT, NT>::fromCOO(coo, [](auto x) { return NT(1); });
}

template <class IT, class NT>
std::array<space::Columnar<IT, NT>, 3> triangleWorstCase(size_t n) {
  std::vector<std::array<IT, 2>> R, S, T;

  for (size_t i = 0; i < n; i++) {
    // N^2 intermediate
    R.push_back({i, 1});
    S.push_back({1, i});
    // But only N triangles
    T.push_back({i, i});
  }

  return {space::Columnar<IT, NT>::fromCOO(R, [](auto x) { return NT(1); }),
          space::Columnar<IT, NT>::fromCOO(S, [](auto x) { return NT(1); }),
          space::Columnar<IT, NT>::fromCOO(T, [](auto x) { return NT(1); })};
}
