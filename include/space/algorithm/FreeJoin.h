#pragma once

#include <array>

#include <space/format/Trie.h>
#include <space/util/Fn.h>

namespace space {

namespace jit {

using namespace fn;

using Var = size_t;
using Rel = size_t;

template <typename... Seqs> struct Nodes;

template <Rel... _Head, typename... _Tail>
struct Nodes<List<_Head...>, _Tail...> {
  using Head = List<_Head...>;
  using Tail = Nodes<_Tail...>;
};

template <Rel... _Head> struct Nodes<List<_Head...>> {
  using Head = List<_Head...>;
  using Tail = Nil;
};

template <class IT, class NT, typename Variables, typename Nodes,
          size_t Batch = 1>
struct FreeJoin;

// Base case: no more nodes
template <class IT, class NT, Var... Vs, size_t Batch>
struct FreeJoin<IT, NT, List<Vs...>, Nil, Batch> {
  static constexpr size_t num_vars = sizeof...(Vs);

  template <class Output, typename Tries>
  static void run(const Tries &tries, std::array<IT, num_vars> &tuple,
                  Value<NT> &value, Output &&output) {
    output(tuple, value);
  }
};

#define FIND_MINIMAL_COVER

// TODO(@altanh): figure out why this is so slow
// Recursive case - batching
template <class IT, class NT, Var... Vs, typename... Ns, size_t Batch>
struct FreeJoin<IT, NT, List<Vs...>, Nodes<Ns...>, Batch> {
  static constexpr size_t num_vars = sizeof...(Vs);
  static constexpr bool last_node = sizeof...(Ns) == 1;

  using Node = typename Nodes<Ns...>::Head;
  template <typename T> using unwrap = std::remove_pointer_t<T>;

  using Tuple = std::array<IT, num_vars>;

  template <class Output, typename Tries>
  static void run(const Tries &all_tries, Tuple &tuple, Value<NT> &value,
                  Output &&output) {
    using std::get;
    using std::vector;

    auto node_tries = fn::Select(all_tries, Node{});

#ifdef FIND_MINIMAL_COVER
    // find minimal cover
    constexpr size_t cover_vars = unwrap<
        std::tuple_element_t<0, decltype(node_tries)>>::Variables::size();
    size_t min_cover_size = std::numeric_limits<size_t>::max();
    size_t min_cover_index = 0;
    fn::ApplyEnumerate(node_tries, [&](auto i, auto trie) {
      // only care if it's a cover
      if constexpr (unwrap<decltype(trie)>::Variables::size() == cover_vars) {
        const size_t cover_size = trie->estimatedSize();
        if (cover_size < min_cover_size) {
          min_cover_size = cover_size;
          min_cover_index = i;
        }
      }
    });
#else
    size_t min_cover_index = 0;
#endif

    auto subtries_init = fn::Apply(node_tries, [](auto trie) {
      return static_cast<typename unwrap<decltype(trie)>::ChildTrie *>(nullptr);
    });

    auto new_tries = fn::ApplyEnumerate(all_tries, [&](auto i, auto trie) {
      if constexpr (fn::Contains<i, Node>::v) {
        return static_cast<typename unwrap<decltype(trie)>::ChildTrie *>(
            nullptr);
      } else {
        return trie;
      }
    });

    // using BatchItem = std::tuple<Tuple, Value<NT>, decltype(subtries_init)>;
    using BatchItem = std::tuple<Tuple, Value<NT>>;
    constexpr size_t TUPLE = 0;
    constexpr size_t VALUE = 1;
    // constexpr size_t SUBTRIES = 2;

    // FIXME: could blow the stack if Batch is too large, but vector is dog slow
    std::array<BatchItem, Batch> batch_items;
    std::array<decltype(subtries_init), Batch> batch_subtries;

    // I wish there was a better way to do this...
    fn::ApplyAtDyn(node_tries, min_cover_index, [&](auto ci, auto ctrie) {
      // iterate over cover
      const size_t iter_size = ctrie->iterSize();
      // iterate in batches
      for (size_t batch = 0; batch < iter_size; batch += Batch) {
        size_t num_items = std::min(Batch, iter_size - batch);
        // populate batch items (tuple, value, new_tries)
        for (size_t b = 0; b < num_items; ++b) {
          typename unwrap<decltype(ctrie)>::Tuple cover_tup;
          auto &batch_item = batch_items[b];
          auto child =
              ctrie->iter(batch + b, cover_tup, get<VALUE>(batch_item));
          fn::Scatter(get<TUPLE>(batch_item), cover_tup,
                      typename unwrap<decltype(ctrie)>::Variables());
          if constexpr (unwrap<decltype(ctrie)>::lastLevel()) {
            if constexpr (!std::is_void_v<NT>) {
              get<VALUE>(batch_item).some *= value.some;
            }
          }
          if constexpr (!last_node) {
            // get<ci>(get<SUBTRIES>(batch_item)) = child;
            get<ci>(batch_subtries[b]) = child;
          }
        }
        // probe the rest
        fn::ApplyEnumerateUntil(node_tries, [&](auto i, auto trie) {
          if constexpr (i != ci) {
            // while loop since we'll be modifying num_items (removing
            // misses)
            size_t b = 0;
            while (b < num_items) {
              typename unwrap<decltype(trie)>::Tuple key;
              auto &batch_item = batch_items[b];
              fn::Gather(key, get<TUPLE>(batch_item),
                         typename unwrap<decltype(trie)>::Variables());
              auto subtrie = trie->lookup(key);
              if (!subtrie) {
                // miss - swap with last element
                std::swap(batch_items[b], batch_items[num_items - 1]);
                --num_items;
              } else {
                if constexpr (unwrap<decltype(trie)>::lastLevel()) {
                  if constexpr (!std::is_void_v<NT>) {
                    get<VALUE>(batch_item).some *= subtrie->value().some;
                  }
                }
                if constexpr (!last_node) {
                  // get<i>(get<SUBTRIES>(batch_item)) = subtrie;
                  get<i>(batch_subtries[b]) = subtrie;
                }
                ++b;
              }
            }
            return num_items == 0;
          }
          return false; // skip cover
        });

        // recurse over batch
        for (size_t b = 0; b < num_items; ++b) {
          auto &batch_item = batch_items[b];
          if constexpr (!last_node) {
            fn::Copy(new_tries, batch_subtries[b], Node{},
                     CountTo<Node::size()>{});
          }
          // recurse
          FreeJoin<IT, NT, List<Vs...>, typename Nodes<Ns...>::Tail,
                   Batch>::run(new_tries, get<TUPLE>(batch_item),
                               get<VALUE>(batch_item),
                               std::forward<Output>(output));
        }
      }
    });
  }
};

// Recursive case - non-batching
template <class IT, class NT, Var... Vs, typename... Ns>
struct FreeJoin<IT, NT, List<Vs...>, Nodes<Ns...>, 1> {
  static constexpr size_t num_vars = sizeof...(Vs);

  using Node = typename Nodes<Ns...>::Head;
  template <typename T> using unwrap = std::remove_pointer_t<T>;

  template <class Output, typename Tries>
  static void run(const Tries &all_tries, std::array<IT, num_vars> &tuple,
                  const Value<NT> &value, Output &&output) {
    auto node_tries = fn::Select(all_tries, Node{});

#ifdef FIND_MINIMAL_COVER
    // find minimal cover
    constexpr size_t cover_vars = unwrap<
        std::tuple_element_t<0, decltype(node_tries)>>::Variables::size();
    size_t min_cover_size = std::numeric_limits<size_t>::max();
    size_t min_cover_index = 0;
    fn::ApplyEnumerate(node_tries, [&](auto i, auto trie) {
      // only care if it's a cover
      if constexpr (unwrap<decltype(trie)>::Variables::size() == cover_vars) {
        const size_t cover_size = trie->estimatedSize();
        if (cover_size < min_cover_size) {
          min_cover_size = cover_size;
          min_cover_index = i;
        }
      }
    });
#else
    size_t min_cover_index = 0;
#endif

    auto new_tries = fn::ApplyEnumerate(all_tries, [&](auto i, auto trie) {
      if constexpr (fn::Contains<i, Node>::v) {
        return static_cast<typename unwrap<decltype(trie)>::ChildTrie *>(
            nullptr);
      } else {
        return trie;
      }
    });

    auto subtries = fn::Apply(node_tries, [](auto trie) {
      return static_cast<typename unwrap<decltype(trie)>::ChildTrie *>(nullptr);
    });

    // I wish there was a better way to do this...
    fn::ApplyAtDyn(node_tries, min_cover_index, [&](auto ci, auto ctrie) {
      // iterate over cover
      const size_t iter_size = ctrie->iterSize();
      for (size_t item = 0; item < iter_size; ++item) {
        Value<NT> v;
        typename unwrap<decltype(ctrie)>::Tuple cover_tup;
        auto child = ctrie->iter(item, cover_tup, v);
        fn::Scatter(tuple, cover_tup,
                    typename unwrap<decltype(ctrie)>::Variables());
        if (ctrie->lastLevel()) {
          if constexpr (!std::is_void_v<NT>) {
            v.some *= value.some;
          }
        }
        std::get<ci>(subtries) = child;
        // probe the rest
        if (fn::ApplyEnumerateUntil(node_tries, [&](auto i, auto trie) {
              if constexpr (i != ci) {
                typename unwrap<decltype(trie)>::Tuple key;
                fn::Gather(key, tuple,
                           typename unwrap<decltype(trie)>::Variables());
                auto subtrie = trie->lookup(key);
                if (!subtrie) {
                  return true; // miss
                }
                if constexpr (unwrap<decltype(trie)>::lastLevel()) {
                  if constexpr (!std::is_void_v<NT>) {
                    v.some *= subtrie->value().some;
                  }
                }
                std::get<i>(subtries) = subtrie;
                return false;
              } else {
                return false; // skip cover
              }
            })) {
          // miss
          continue;
        }
        // merge subtries
        fn::Copy(new_tries, subtries, Node{}, CountTo<Node::size()>{});
        // recurse
        FreeJoin<IT, NT, List<Vs...>, typename Nodes<Ns...>::Tail>::run(
            new_tries, tuple, v, std::forward<Output>(output));
      }
    });
  }
};

} // namespace jit

} // namespace space
