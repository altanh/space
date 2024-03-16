#pragma once

#include <algorithm>
#include <cstdint>
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
  Columnar<IT, NT> *data; // TODO: shared_ptr, probably
};

// TODO: put this somewhere else?
template <class NTOut> struct Value;
template <> struct Value<void> {
  template <class NT> static Value<void> make(NT init) { return {}; }
};
template <class NTOut> struct Value {
  NTOut some;

  Value() : some(NTOut(1)) {}

  template <class NT> static Value<NTOut> make(NT init) {
    return {static_cast<NTOut>(init)};
  }
};

using Dim = uint8_t;

template <class IT, class NT> class LazyTrie {
  // cover is true if the last level of the trie is a cover for its node in the
  // join plan
  bool is_cover;
  // is_cover && nlevels == 1
  bool is_last_cover;
  // number of levels in the trie; 0 means the trie is a leaf, i.e. nvars = 0
  Dim nlevels;
  // nvars[0]; see below for nvars
  Dim this_nvars;
  // nvars[i] is the number of variables at level i
  Dim *nvars;
  // perm[i] is conceptually an array of size nvars[i] that maps the trie
  // variables (by index) to the relation variables (by index). For example,
  // given R(x, y, z) and a trie with schema [[y, z], [x]], perm = [[1, 2],
  // [0]]. In reality, perm is stored as a flat array, so perm[0] = 1, perm[1] =
  // 2, perm[2] = 0.
  Dim *perm;
  Relation<IT, NT> *relation;
  // to avoid an unnecessary indirection
  Columnar<IT, NT> *data; // I should probably just use Rust

  // TODO(@altanh): experiment with storing the Trie directly in the HashTable,
  //                can avoid some indirections. However, this could make
  //                resizing more costly.
  // std::unique_ptr<HashTable<IT, std::unique_ptr<LazyTrie<IT, NT>>>> map;
  std::unique_ptr<HashTable<IT, LazyTrie<IT, NT>>> map;
  std::unique_ptr<std::vector<size_t>> vec;

public:
  LazyTrie() = default;

  LazyTrie(bool root, bool is_cover, Dim nlevels, Dim *nvars, Dim *perm,
           Relation<IT, NT> *relation)
      : is_cover(is_cover), is_last_cover(is_cover && nlevels == 1),
        nlevels(nlevels), this_nvars(nvars[0]), nvars(nvars), perm(perm),
        relation(relation), data(relation->data) {
    if (root) {
      vec = nullptr; // root level can avoid materializing [0, ..., N-1]
    } else {
      vec = std::make_unique<std::vector<size_t>>();
    }
  }

  static LazyTrie<IT, NT> fromSchema(Relation<IT, NT> *rel,
                                     const Schema &schema, bool is_cover,
                                     std::unique_ptr<Dim[]> *nvars_out,
                                     std::unique_ptr<Dim[]> *perm_out) {
    // build var -> idx map
    std::unordered_map<Variable, Dim> var_to_idx;
    for (Dim i = 0; i < rel->vars.size(); i++) {
      var_to_idx[rel->vars[i]] = i;
    }
    // calculate nvars
    Dim nlevels = schema.size();
    std::unique_ptr<Dim[]> nvars = std::make_unique<Dim[]>(nlevels + 1);
    Dim nvars_total = 0;
    for (Dim i = 0; i < nlevels; i++) {
      nvars_total += schema[i].size();
      nvars[i] = schema[i].size();
    }
    nvars[nlevels] = 0; // last level is a leaf
    // build permutation map
    std::unique_ptr<Dim[]> perm = std::make_unique<Dim[]>(nvars_total);
    Dim offset = 0;
    for (Dim level = 0; level < nlevels; level++) {
      for (Dim var = 0; var < schema[level].size(); var++) {
        perm[offset + var] = var_to_idx[schema[level][var]];
      }
      offset += nvars[level];
    }
    *nvars_out = std::move(nvars);
    *perm_out = std::move(perm);
    return LazyTrie<IT, NT>(true, is_cover, nlevels, nvars_out->get(),
                            perm_out->get(), rel);
  }

  template <class ITOut, class NTOut>
  void iter(size_t idx, ITOut *out, Value<NTOut> *val) {
    // TODO: will we ever need to fetch the value for a map?
    if (map) {
      map->iter(idx, out);
    } else if (is_last_cover) {
      // last subatom, and a cover of the node
      size_t item = vec ? (*vec)[idx] : idx;
      for (Dim var = 0; var < this_nvars; var++) {
        out[var] = data->dim(perm[var])[item];
      }
      if constexpr (!(std::is_void<NTOut>::value || std::is_void<NT>::value)) {
        val->some = data->val(item)[0];
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
    // return result ? result->get() : nullptr;
    return result;
  }

  // WARNING: this will force if necessary!
  size_t getIterSize() {
    if (map) {
      return map->getSize();
    } else if (is_last_cover) {
      return vec ? vec->size() : data->getSize();
    } else {
      force();
      return map->getSize();
    }
  }

  // WARNING: this will force if necessary!
  size_t getMapSize() {
    if (!map) {
      force();
    }
    return map->getSize();
  }

  // TODO: improvements to this?
  size_t getEstimatedSize() const {
    if (map) {
      return map->getSize();
    } else {
      return vec ? vec->size() : data->getSize();
    }
  }

  size_t getNumVars() const { return this_nvars; }

  const Relation<IT, NT> &getRelation() const { return *relation; }

  bool isLastCover() const { return is_last_cover; }

  Dim getNumLevels() const { return nlevels; }

  void printSchema() const {
    std::cout << relation->name << ": ";
    Dim lv = nlevels;
    Dim *nv = nvars;
    Dim *p = perm;
    while (true) {
      std::cout << "[";
      for (Dim i = 0; i < *nv; i++) {
        std::cout << relation->vars[*p];
        if (i < *nv - 1) {
          std::cout << ", ";
        }
        p++;
      }
      std::cout << "]";
      if (lv == 0) {
        break;
      }
      std::cout << ", ";
      lv--;
      nv++;
    }
    std::cout << std::endl;
  }

protected:
  void force() {
    if (map) {
      return;
    }
    // allocate some stack space for the (small) keys
    IT *key = static_cast<IT *>(alloca(this_nvars * sizeof(IT)));

    // map = std::make_unique<HashTable<IT, std::unique_ptr<LazyTrie<IT, NT>>>>(
    //     this_nvars);
    map = std::make_unique<HashTable<IT, LazyTrie<IT, NT>>>(this_nvars);

    // root case
    if (!vec) {
      for (size_t index = 0; index < data->getSize(); index++) {
        for (Dim var = 0; var < this_nvars; var++) {
          key[var] = data->dim(perm[var])[index];
        }
        if (map->lookup(key) == nullptr) {
          // map->emplace(key, std::make_unique<LazyTrie<IT, NT>>(
          //                       false, is_cover, nlevels - 1, nvars + 1,
          //                       perm + this_nvars, relation));
          map->emplace(key,
                       LazyTrie<IT, NT>(false, is_cover, nlevels - 1, nvars + 1,
                                        perm + this_nvars, relation));
        }
        // (*map->lookup(key))->vec->push_back(index);
        map->lookup(key)->vec->push_back(index);
      }
    } else {
      for (size_t index : *vec) {
        for (Dim var = 0; var < this_nvars; var++) {
          key[var] = data->dim(perm[var])[index];
        }
        if (map->lookup(key) == nullptr) {
          // map->emplace(key, std::make_unique<LazyTrie<IT, NT>>(
          //                       false, is_cover, nlevels - 1, nvars + 1,
          //                       perm + this_nvars, relation));
          map->emplace(key,
                       LazyTrie<IT, NT>(false, is_cover, nlevels - 1, nvars + 1,
                                        perm + this_nvars, relation));
        }
        // (*map->lookup(key))->vec->push_back(index);
        map->lookup(key)->vec->push_back(index);
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
