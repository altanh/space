#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "Columnar.h"
#include "HashTable.h"

namespace space {

// COLT from Free Join [Wang et al. 2023]
using Variable = std::string;
using Schema = std::vector<std::vector<Variable>>;

template <class IT, class NT> struct Relation {
  Variable name;
  std::vector<Variable> vars;
  std::shared_ptr<Columnar<IT, NT>> data; // TODO: shared_ptr, probably
};

// TODO: put this somewhere else?
template <class NTOut> struct Value;
template <> struct Value<void> {};
template <class NTOut> struct Value {
  NTOut some;
};

template <class IT, class NT> class LazyTrie {
  Relation<IT, NT> relation;
  bool is_suffix;
  Schema schema;
  std::vector<Variable> vars;
  size_t nvars;

  std::unique_ptr<size_t[]> dim_keys;

  // TODO(@altanh): experiment with storing the Trie directly in the HashTable,
  //                can avoid some indirections
  std::unique_ptr<HashTable<IT, std::unique_ptr<LazyTrie<IT, NT>>>> map;
  std::unique_ptr<std::vector<size_t>> vec;

public:
  // Top-level Trie case
  LazyTrie(Relation<IT, NT> rel, Schema schema) : LazyTrie(rel, schema, true) {}

  LazyTrie(Relation<IT, NT> rel, Schema schema, bool root) {
    std::cout << "Making Trie with schema: ";
    for (const auto &vars : schema) {
      std::cout << "[";
      if (vars.size() > 0) {
        std::cout << vars[0];
      }
      for (size_t i = 1; i < vars.size(); ++i) {
        std::cout << ", " << vars[i];
      }
      std::cout << "]";
    }
    std::cout << std::endl;
    this->relation = rel;
    this->schema = schema;
    this->vars = schema.empty() ? std::vector<Variable>{} : schema[0];
    this->nvars = vars.size();
    this->is_suffix = _isSuffix(vars, rel.vars);

    map = nullptr;
    if (root) {
      vec = nullptr; // root level can avoid materializing [0, ..., N-1]
    } else {
      vec = std::make_unique<std::vector<size_t>>();
    }

    dim_keys = std::make_unique<size_t[]>(nvars);
    if (nvars > 0) {
      size_t var = 0;
      for (size_t i = 0; i < relation.vars.size(); i++) {
        if (relation.vars[i] == vars[var]) {
          dim_keys[var] = i;
          var++;
          if (var == nvars) {
            break;
          }
        }
      }
    }
  }

  template <class ITOut, class NTOut>
  void iter(size_t idx, ITOut *out, Value<NTOut> *val) {
    // TODO: will we ever need to fetch the value for a map?
    if (map) {
      map->iter(idx, out);
    } else if (is_suffix) {
      size_t item = vec ? (*vec)[idx] : idx;
      for (size_t i = 0; i < nvars; i++) {
        out[i] = relation.data->dim(dim_keys[i])[item];
      }
      if constexpr (!(std::is_void<NTOut>::value || std::is_void<NT>::value)) {
        val->some = *relation.data->val(item);
      }
    } else {
      force();
      map->iter(idx, out);
    }
  }

  LazyTrie<IT, NT> *lookup(IT *key) {
    if (!map) {
      force();
    }
    auto result = map->lookup(key);
    return result ? result->get() : nullptr;
  }

  // WARNING: this will force if necessary!
  size_t getSize() {
    if (map) {
      return map->getSize();
    } else if (is_suffix) {
      return vec ? vec->size() : relation.data->getSize();
    } else {
      force();
      return map->getSize();
    }
  }

  size_t getNumVars() const { return nvars; }

  const Schema &getSchema() const { return schema; }

  const Relation<IT, NT> &getRelation() const { return relation; }

  bool isSuffix() const { return is_suffix; }

protected:
  static bool _isSuffix(std::vector<Variable> smaller,
                        std::vector<Variable> larger) {
    const size_t n = smaller.size();
    const size_t m = larger.size();
    for (size_t i = 0; i < n; i++) {
      if (smaller[n - i - 1] != larger[m - i - 1]) {
        return false;
      }
    }
    return true;
  }

  void force() {
    if (map) {
      return;
    }
    // allocate some stack space for the (small) keys
    IT *key = static_cast<IT *>(alloca(nvars * sizeof(IT)));

    // map = new HashTable<IT, LazyTrie<IT, NT> *>(nvars);
    map = std::make_unique<HashTable<IT, std::unique_ptr<LazyTrie<IT, NT>>>>(
        nvars);
    // root case
    if (!vec) {
      for (size_t i = 0; i < relation.data->getSize(); i++) {
        for (size_t j = 0; j < nvars; j++) {
          key[j] = relation.data->dim(dim_keys[j])[i];
        }
        if (map->lookup(key) == nullptr) {
          // TODO: make this faster by not copying the schema
          Schema new_schema(schema.begin() + 1, schema.end());
          map->emplace(key, std::make_unique<LazyTrie<IT, NT>>(
                                relation, new_schema, false));
        }
        (*map->lookup(key))->vec->push_back(i);
      }
    } else {
      for (size_t i : *vec) {
        for (size_t j = 0; j < nvars; j++) {
          key[j] = relation.data->dim(dim_keys[j])[i];
        }
        if (map->lookup(key) == nullptr) {
          // TODO
          Schema new_schema(schema.begin() + 1, schema.end());
          map->emplace(key, std::make_unique<LazyTrie<IT, NT>>(
                                relation, new_schema, false));
        }
        (*map->lookup(key))->vec->push_back(i);
      }
      vec = nullptr;
    }
  }
};

// Type-erased list of tries parametrized by output types
// template <class ITOut, class NTOut> class TrieList {
//   // Interface
//   class TrieInterface {
//   public:
//     virtual void iter(size_t idx, ITOut *out, Value<NTOut> *val) = 0;
//     virtual TrieInterface *lookup(ITOut *key) = 0;
//     virtual size_t getSize() = 0;
//     virtual size_t getNumVars() const = 0;
//     virtual const Schema &getSchema() const = 0;
//     virtual bool isSuffix() const = 0;
//   };

//   // Derived type
//   template <class IT, class NT> class TrieImpl : public TrieInterface {
//     LazyTrie<IT, NT> *trie;

//   public:
//     TrieImpl(LazyTrie<IT, NT> *trie) : trie(trie) {}

//     void iter(size_t idx, ITOut *out, Value<NTOut> *val) override {
//       trie->iter(idx, out, val);
//     }

//     TrieInterface *lookup(ITOut *key) override {
//       return new TrieImpl<IT, NT>(trie->lookup(key));
//     }
//   }
// };

// A(i,j) stored as CSR:
// Trie(A, [[i], [j]])
//
// A(i,j) stored as CSC:
// Trie(A, [[j], [i]])

// Data lifecycle:
// 1. Tensor data starts in columnar storage (sort of "transposed" COO)
// 2. COO is read into Trie according to query plan

} // namespace space
