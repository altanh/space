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

template <class IT, class NT, typename Variables, typename Nodes>
struct FreeJoin;

// Base case: no more nodes
template <class IT, class NT, size_t... Vs>
struct FreeJoin<IT, NT, List<Vs...>, Nil> {
  static constexpr size_t num_vars = sizeof...(Vs);

  template <class Output, typename Tries>
  static void run(const Tries &tries, std::array<IT, num_vars> &tuple,
                  Value<NT> &value, Output &&output) {
    output(tuple, value);
  }
};

// Recursive case
template <class IT, class NT, size_t... Vs, typename... Ns>
struct FreeJoin<IT, NT, List<Vs...>, Nodes<Ns...>> {
  static constexpr size_t num_vars = sizeof...(Vs);

  using Node = typename Nodes<Ns...>::Head;
  template <typename T> using unwrap = std::remove_pointer_t<T>;

  template <class Output, typename Tries>
  static void run(const Tries &all_tries, std::array<IT, num_vars> &tuple,
                  const Value<NT> &value, Output &&output) {
    auto node_tries = fn::Select(all_tries, Node{});

#define FIND_MINIMAL_COVER

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
        ctrie->iter(item, cover_tup, v);
        fn::Scatter(tuple, cover_tup,
                    typename unwrap<decltype(ctrie)>::Variables());
        if (ctrie->lastLevel()) {
          if constexpr (!std::is_void_v<NT>) {
            v.some *= value.some;
          }
          std::get<ci>(subtries) = nullptr;
        } else {
          std::get<ci>(subtries) = ctrie->lookup(cover_tup);
        }
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
