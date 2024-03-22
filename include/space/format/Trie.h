#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "Columnar.h"
#include "HashTable.h"
#include <space/util/Fn.h>

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

namespace jit {

using namespace fn;

template <class IT, class NT, typename Schema, typename Variables, size_t Rank,
          bool Root = true, size_t Bound = 0>
class LazyTrie;

// Base case: empty schema
template <class IT, class NT, bool Root, size_t Rank, size_t Bound,
          size_t... Vs>
class LazyTrie<IT, NT, Nil, List<Vs...>, Rank, Root, Bound> {
public:
  using Variables = Nil;
  using Tuple = std::array<IT, 0>;

private:
  template <class IT2, class NT2, typename S, typename V, size_t Ra, bool Ro,
            size_t B>
  friend class LazyTrie;

  Value<NT> value_;

public:
  LazyTrie() = default;
  LazyTrie(LazyTrie &&other) = default;
  LazyTrie &operator=(LazyTrie &&other) = default;
  LazyTrie(Columnar<IT, NT, Rank> *data) = delete;
  LazyTrie(Columnar<IT, NT, Rank> *data, size_t index) {
    if constexpr (!std::is_void_v<NT>) {
      value_.some = data->value(index);
    }
  }

  void push(size_t index) {
    // should be unreachable
    throw std::runtime_error("push called on leaf trie");
  }

  template <class ITOut, class NTOut>
  void iter(size_t idx, Tuple &out, Value<NTOut> &val) {
    // idx unused, only one element in leaf tries
    if constexpr (!(std::is_void_v<NTOut> || std::is_void_v<NT>)) {
      val.some = static_cast<NTOut>(value_.some);
    }
  }

  Value<NT> &value() { return value_; }

  size_t iterSize() { return 1; }

  size_t mapSize() { return 0; }

  size_t estimatedSize() { return 1; }

  void print(size_t level = 0) {
    for (size_t j = 0; j < level; ++j) {
      std::cout << "    ";
    }
    std::cout << "()";
    if constexpr (!std::is_void_v<NT>) {
      std::cout << " -> " << value_.some;
    }
    std::cout << std::endl;
  }

  void materialize() {}
};

// Recursive case

// template <class IT, class NT, typename Schema, typename Variables, size_t
// Rank,
//           bool Root = true>
// class LazyTrie;

template <class IT, class NT, size_t Rank, size_t Bound, bool Root, size_t N,
          size_t... Ns, size_t... Vs>
class LazyTrie<IT, NT, List<N, Ns...>, List<Vs...>, Rank, Root, Bound> {
public:
  using Variables = typename Take<N, List<Vs...>>::T;
  using Tuple = std::array<IT, N>;
  static_assert(Variables::size() == N, "Variables size mismatch");
  using ChildSchema = List<Ns...>;
  using ChildVariables = typename Tail<N, List<Vs...>>::T;
  using ChildTrie =
      LazyTrie<IT, NT, ChildSchema, ChildVariables, Rank, false, Bound + N>;

private:
  template <class IT2, class NT2, typename S, typename V, size_t Ra, bool Ro,
            size_t Re>
  friend class LazyTrie;

  using ProjectionSeq = typename Range<Bound, Bound + N>::T;

  using Map = StaticHashTable<IT, ChildTrie, N>;

  // Could use a variant, but might honestly be clunkier.
  std::unique_ptr<Map> map_;
  std::unique_ptr<std::vector<size_t>> vec_;

  Columnar<IT, NT, Rank> *data_;

  template <size_t... Is, size_t... Cs>
  void writeTuple(size_t item, std::array<IT, N> &out, List<Is...>,
                  List<Cs...>) {
    (..., (out[Is] = data_->column(Cs)[item]));
  }

  void writeTuple(size_t item, std::array<IT, N> &out) {
    writeTuple(item, out, CountTo<N>{}, ProjectionSeq{});
  }

public:
  LazyTrie() = default;

  LazyTrie(LazyTrie &&other) = default;

  // Move assignment
  LazyTrie &operator=(LazyTrie &&other) = default;

  LazyTrie(Columnar<IT, NT, Rank> *data) : data_(data) {
    map_ = nullptr;
    if constexpr (Root) {
      vec_ = nullptr;
    } else {
      vec_ = std::make_unique<std::vector<size_t>>();
    }
  }

  LazyTrie(Columnar<IT, NT, Rank> *data, size_t index) : data_(data) {
    vec_ = std::make_unique<std::vector<size_t>>();
    vec_->push_back(index);
  }

  template <class ITOut, class NTOut>
  void iter(size_t idx, std::array<ITOut, N> &out, Value<NTOut> &val) {
    if (map_) {
      map_->iter(idx, out);
    } else if constexpr (lastLevel()) {
      size_t item = Root ? idx : (*vec_)[idx];
      writeTuple(item, out);
      if constexpr (!(std::is_void_v<NTOut> || std::is_void_v<NT>)) {
        val.some = data_->value(item);
      }
    } else {
      force();
      map_->iter(idx, out);
    }
  }

  ChildTrie *lookup(const Tuple &key) {
    if (!map_) {
      force();
    }
    return map_->lookup(key);
  }

  size_t iterSize() {
    if (map_) {
      return map_->getSize();
    } else if constexpr (lastLevel()) {
      return Root ? data_->size() : vec_->size();
    } else {
      force();
      return map_->getSize();
    }
  }

  size_t mapSize() {
    if (!map_) {
      force();
    }
    return map_->getSize();
  }

  size_t estimatedSize() const {
    if (map_) {
      return map_->getSize();
    } else {
      return Root ? data_->size() : vec_->size();
    }
  }

  void print(size_t level = 0) {
    const size_t size = iterSize();
    Tuple key;
    for (size_t i = 0; i < size; ++i) {
      Value<NT> value;
      iter(i, key, value);
      for (size_t j = 0; j < level; ++j) {
        std::cout << "    ";
      }
      std::cout << "(";
      if (N > 0) {
        std::cout << key[0];
      }
      for (size_t j = 1; j < N; ++j) {
        std::cout << ", " << key[j];
      }
      std::cout << ")";
      if constexpr (lastLevel()) {
        std::cout << " -> " << value.some << std::endl;
      } else {
        std::cout << std::endl;
        lookup(key)->print(level + 1);
      }
    }
  }

  static constexpr bool lastLevel() { return ChildSchema::size() == 0; }

  void materialize() {
    if (map_) {
      return;
    }

    Tuple key;
    map_ = std::make_unique<Map>();
    if (Root) {
      for (size_t index = 0; index < data_->size(); index++) {
        writeTuple(index, key);
        ChildTrie *child = map_->lookup(key);
        if (child) {
          child->push(index);
        } else {
          map_->emplace(key, ChildTrie(data_, index));
        }
      }
    } else {
      for (size_t index : *vec_) {
        writeTuple(index, key);
        ChildTrie *child = map_->lookup(key);
        if (child) {
          child->push(index);
        } else {
          map_->emplace(key, ChildTrie(data_, index));
        }
      }
      vec_ = nullptr;
      // iterate over children and materialize
      const size_t map_size = map_->getSize();
      for (size_t i = 0; i < map_size; i++) {
        map_->iter(i)->materialize();
      }
    }
  }

protected:
  void push(size_t index) { vec_->push_back(index); }

  void force() {
    if (map_) {
      return;
    }

    Tuple key;
    map_ = std::make_unique<Map>();
    if (Root) {
      for (size_t index = 0; index < data_->size(); index++) {
        writeTuple(index, key);
        ChildTrie *child = map_->lookup(key);
        if (child) {
          child->push(index);
        } else {
          map_->emplace(key, ChildTrie(data_, index));
        }
      }
    } else {
      for (size_t index : *vec_) {
        writeTuple(index, key);
        ChildTrie *child = map_->lookup(key);
        if (child) {
          child->push(index);
        } else {
          map_->emplace(key, ChildTrie(data_, index));
        }
      }
      vec_ = nullptr;
    }
  }
};

} // namespace jit

// Trie with compile-time schema, provided as a parameter pack of ints
// Schema[i] is the number of variables at level i.
// Perm[j] has length sum(Schema), and maps trie variable j to the relation.
template <class IT, class NT, bool Root, bool LeafCover, typename Schema,
          typename Perm, typename OutPerm = Perm>
class StaticLazyTrie;

// Base case: no variables
template <class IT, class NT, bool Root, bool LeafCover, size_t... Ps,
          size_t... OPs>
class StaticLazyTrie<IT, NT, Root, LeafCover, std::index_sequence<>,
                     std::index_sequence<Ps...>, std::index_sequence<OPs...>> {
public:
  // Expose template parameters.
  using OutPerm = std::index_sequence<>;
  using Tuple = std::array<IT, 0>;

private:
  // Make friends.
  template <class IT2, class NT2, bool R, bool LC, typename S, typename P,
            typename OP2>
  friend class StaticLazyTrie;

  Columnar<IT, NT> *data;

  template <typename T> using ptr = std::unique_ptr<T>;
  template <typename T> using vector = std::vector<T>;
  template <typename T, size_t M> using array = std::array<T, M>;
  template <size_t... is> using iseq = std::index_sequence<is...>;

  // No map needed, only a single index is stored.
  size_t index;

public:
  StaticLazyTrie() = default;

  StaticLazyTrie(StaticLazyTrie &&other) = default;

  StaticLazyTrie &operator=(StaticLazyTrie &&other) = default;

  StaticLazyTrie(Columnar<IT, NT> *data) : data(data) {}

  template <class ITOut, class NTOut>
  void iter(size_t idx, array<ITOut, 0> &out, Value<NTOut> &val) {
    if constexpr (!(std::is_void<NTOut>::value || std::is_void<NT>::value)) {
      val.some = *data->val(index);
    }
  }

  size_t getIterSize() { return 1; }

  size_t getMapSize() { return 0; }

  size_t getEstimatedSize() const { return 1; }

  void push(size_t index) { this->index = index; }

  void print(size_t level = 0) {
    for (size_t j = 0; j < level; ++j) {
      std::cout << "    ";
    }
    std::cout << "()";
    if constexpr (!(std::is_void<NT>::value)) {
      std::cout << " -> " << *data->val(index);
    }
    std::cout << std::endl;
  }
};

// helper to get the tail of an index_sequence from some element onwards
template <size_t Index, typename Seq> struct tail;

template <size_t Head, size_t... Tail>
struct tail<0, std::index_sequence<Head, Tail...>> {
  using type = std::index_sequence<Head, Tail...>;
};

template <> struct tail<0, std::index_sequence<>> {
  using type = std::index_sequence<>;
};

template <size_t Index, size_t Head, size_t... Tail>
struct tail<Index, std::index_sequence<Head, Tail...>> {
  using type = typename tail<Index - 1, std::index_sequence<Tail...>>::type;
};

// helper to concatenate two index_sequences
template <typename Seq1, typename Seq2> struct concat;

template <size_t... Is1, size_t... Is2>
struct concat<std::index_sequence<Is1...>, std::index_sequence<Is2...>> {
  using type = std::index_sequence<Is1..., Is2...>;
};

// helper to get the first n elements of an index_sequence
template <size_t N, typename Seq> struct take;

template <size_t N, size_t Head, size_t... Tail>
struct take<N, std::index_sequence<Head, Tail...>> {
  using type = typename concat<
      std::index_sequence<Head>,
      typename take<N - 1, std::index_sequence<Tail...>>::type>::type;
};

template <size_t Head, size_t... Tail>
struct take<1, std::index_sequence<Head, Tail...>> {
  using type = std::index_sequence<Head>;
};

template <> struct take<0, std::index_sequence<>> {
  using type = std::index_sequence<>;
};

template <class IT, class NT, bool Root, bool LeafCover, size_t N, size_t... Ns,
          size_t... Ps, size_t... OPs>
class StaticLazyTrie<IT, NT, Root, LeafCover, std::index_sequence<N, Ns...>,
                     std::index_sequence<Ps...>, std::index_sequence<OPs...>> {
public:
  // Expose template parameters.
  static constexpr bool root = Root;
  static constexpr bool leaf_cover = LeafCover;
  static constexpr size_t num_vars = N;
  static constexpr bool last_level = sizeof...(Ns) == 0;

  using OutPerm = typename take<N, std::index_sequence<OPs...>>::type;

  using ChildPerm = typename tail<N, std::index_sequence<Ps...>>::type;
  using ChildOutPerm = typename tail<N, std::index_sequence<OPs...>>::type;
  using ChildTrie =
      StaticLazyTrie<IT, NT, false, LeafCover, std::index_sequence<Ns...>,
                     ChildPerm, ChildOutPerm>;

  using Tuple = std::array<IT, N>;
  using V = Value<NT>;

private:
  // Make friends.
  template <class IT2, class NT2, bool R, bool LC, typename S, typename P,
            typename OP2>
  friend class StaticLazyTrie;

  Columnar<IT, NT> *data;

  template <typename T> using ptr = std::unique_ptr<T>;
  template <typename T> using vector = std::vector<T>;
  template <typename T, size_t M> using array = std::array<T, M>;

  using Map = StaticHashTable<IT, ChildTrie, N>;

  ptr<StaticHashTable<IT, ChildTrie, N>> map;
  ptr<vector<size_t>> vec;

  // helper function to get the j-th element of Ps
  template <size_t j> static constexpr size_t P() {
    return std::get<j>(std::make_tuple(Ps...));
  }

  // Helper function to write out a tuple from the data.
  // Loop over the N variables and write the data to the out array using P.
  template <size_t... is>
  void writeTuple(size_t item, array<IT, N> &out, std::index_sequence<is...>) {
    (..., (out[is] = data->dim(P<is>())[item]));
  }

public:
  StaticLazyTrie() = default;

  StaticLazyTrie(StaticLazyTrie &&other) = default;

  // Move assignment
  StaticLazyTrie &operator=(StaticLazyTrie &&other) = default;

  StaticLazyTrie(Columnar<IT, NT> *data) : data(data) {
    if constexpr (Root) {
      vec = nullptr;
    } else {
      vec = std::make_unique<vector<size_t>>();
    }
  }

  template <class ITOut, class NTOut>
  void iter(size_t idx, array<ITOut, N> &out, Value<NTOut> &val) {
    if (map) {
      map->iter(idx, out);
    } else if constexpr (last_level && LeafCover) {
      size_t item = Root ? idx : (*vec)[idx];
      writeTuple(item, out, std::make_index_sequence<N>{});
      if constexpr (!(std::is_void<NTOut>::value || std::is_void<NT>::value)) {
        val.some = *data->val(item);
      }
    } else {
      force();
      map->iter(idx, out);
    }
  }

  ChildTrie *lookup(Tuple &key) {
    if (!map) {
      force();
    }
    auto result = map->lookup(key);
    return result;
  }

  void push(size_t index) { vec->push_back(index); }

  size_t getIterSize() {
    if (map) {
      return map->getSize();
    } else if constexpr (last_level && LeafCover) {
      return Root ? data->getSize() : vec->size();
    } else {
      force();
      return map->getSize();
    }
  }

  size_t getMapSize() {
    if (!map) {
      force();
    }
    return map->getSize();
  }

  size_t getEstimatedSize() const {
    if (map) {
      return map->getSize();
    } else {
      return Root ? data->getSize() : vec->size();
    }
  }

  void print(size_t level = 0) {
    const size_t size = getIterSize();
    Tuple key;
    for (size_t i = 0; i < size; ++i) {
      Value<NT> value;
      iter(i, key, value);
      for (size_t j = 0; j < level; ++j) {
        std::cout << "    ";
      }
      std::cout << "(";
      if (N > 0) {
        std::cout << key[0];
      }
      for (size_t j = 1; j < N; ++j) {
        std::cout << ", " << key[j];
      }
      std::cout << ")";
      if constexpr (last_level && LeafCover) {
        std::cout << " -> " << value.some << std::endl;
      } else {
        std::cout << std::endl;
        lookup(key)->print(level + 1);
      }
    }
  }

protected:
  void force() {
    if (map) {
      return;
    }

    Tuple key;
    map = std::make_unique<Map>();
    if (Root) {
      for (size_t index = 0; index < data->getSize(); index++) {
        writeTuple(index, key, std::make_index_sequence<N>{});
        if (map->lookup(key) == nullptr) {
          map->emplace(key, ChildTrie(data));
        }
        map->lookup(key)->push(index);
      }
    } else {
      for (size_t index : *vec) {
        writeTuple(index, key, std::make_index_sequence<N>{});
        if (map->lookup(key) == nullptr) {
          map->emplace(key, ChildTrie(data));
        }
        map->lookup(key)->push(index);
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
