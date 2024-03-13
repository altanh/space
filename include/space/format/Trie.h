#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "Columnar.h"
#include "HashTable.h"

namespace space {

// inspired by TACO and the Generalized Hash Trie from Free Join
using Variable = std::string;
using Schema = std::vector<std::vector<Variable>>;

template <class IT, class NT> struct Relation {
  Variable name;
  std::vector<Variable> vars;
  Columnar<IT, NT> data;
};

template <class IT, class NT> class LazyTrie {
  Relation<IT, NT> relation;
  const Schema schema;
  const std::vector<Variable> vars;
  const size_t nvars;

  HashTable<IT, LazyTrie<IT, NT> *> *map;
  std::vector<size_t> *vec;

public:
  // Top-level Trie case
  Trie(Relation rel, Schema schema) : Trie(rel, schema, true) {}

  Trie(Relation rel, Schema schema, bool root)
      : relation(rel), schema(schema), vars(schema[0]), nvars(vars.size()) {
    map = nullptr;
    if (root) {
      vec = nullptr; // root level can avoid materializing [0, ..., N-1]
    } else {
      vec = new std::vector<size_t>();
    }
  }

  // TODO: iter()

  // TODO: get()

protected:
  void force() {
    if (map) {
      return;
    }

    // allocate some stack space for the (small) keys
    IT *idxs = alloca(nvars * sizeof(size_t));
    IT *key = alloca(nvars * sizeof(IT));

    // get index of trie vars in relation vars for column lookup
    for (size_t i = 0; i < nvars; i++) {
      auto it = std::find(relation.vars.begin(), relation.vars.end(), vars[i]);
      if (it == relation.vars.end()) {
        throw std::runtime_error("variable not found");
      }
      idxs[i] = std::distance(relation.vars.begin(), it);
    }

    map = new HashTable<IT, LazyTrie<IT, NT> *>(nvars);
    // root case
    if (vec == nullptr) {
      for (size_t i = 0; i < relation.data.nnz; i++) {
        for (size_t j = 0; j < nvars; j++) {
          key[j] = *relation.data.dim(idxs[j]);
        }
        if (map->lookup(key) == nullptr) {
          map->insert(key, new LazyTrie<IT, NT>(relation, schema, false);
        }
        map->lookup(key)->vec->push_back(i);
      }
    } else {
      for (size_t i : *vec) {
        for (size_t j = 0; j < nvars; j++) {
          key[j] = *relation.data.dim(idxs[j]);
        }
        if (map->lookup(key) == nullptr) {
          map->insert(key, new LazyTrie<IT, NT>(relation, schema, false);
        }
        map->lookup(key)->vec->push_back(i);
      }
      delete vec;
      vec = nullptr;
    }
  }
};

// A(i,j) stored as CSR:
// Trie(A, [[i], [j]])
//
// A(i,j) stored as CSC:
// Trie(A, [[j], [i]])

// Data lifecycle:
// 1. Tensor data starts in columnar storage (sort of "transposed" COO)
// 2. COO is read into Trie according to query plan

} // namespace space