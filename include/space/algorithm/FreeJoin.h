#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <space/format/Trie.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace space {

namespace jit {

using Id = size_t;

// Free Join inputs:
// - VO: number of output variables
// - VP: number of bound variables in the output
// - Tuple of StaticLazyTrie pointers corresponding to the relations
// - Array of size VO for storing the output tuple
// - List of nodes. A node is just a list of indices into the relations tuple.
//   This list represented at compile time using a list of index_sequences

// NodeList type. Stores a list of index_sequences.
template <typename... Seqs> struct NodeList;

template <Id... H, typename... T>
struct NodeList<std::index_sequence<H...>, T...> {
  using Head = std::index_sequence<H...>;
  using Tail = NodeList<T...>;
};

template <Id... H> struct NodeList<std::index_sequence<H...>> {
  using Head = std::index_sequence<H...>;
  using Tail = NodeList<>;
};

using std::array;
using std::get;
using std::index_sequence;
using std::tuple;

// template metaprogramming stuff
namespace fn {

// helper to apply a function to each trie in the given index sequence
template <typename Seq, typename Ts> struct Apply;

template <Id... Is, typename... Ts>
struct Apply<index_sequence<Is...>, tuple<Ts...>> {
  template <typename F> static void run(const tuple<Ts...> &tries, F &&f) {
    (f(get<Is>(tries), std::integral_constant<Id, Is>()), ...);
  }
};

template <typename Seq, typename Ts> struct ApplyUntil;

template <Id... Is, typename... Ts>
struct ApplyUntil<index_sequence<Is...>, tuple<Ts...>> {
  template <typename F> static bool run(const tuple<Ts...> &tries, F &&f) {
    return (f(get<Is>(tries), std::integral_constant<Id, Is>()) || ...);
  }
};

// helper to map a function over each trie in the given index sequence
template <typename Seq, typename Ts> struct Map;

template <Id... Is, typename... Ts>
struct Map<index_sequence<Is...>, tuple<Ts...>> {
  template <typename F> static auto run(const tuple<Ts...> &tries, F &&f) {
    return std::make_tuple(
        f(get<Is>(tries), std::integral_constant<Id, Is>())...);
  }
};

// helper to check if a number shows up in an index sequence
template <Id I, typename Seq> struct Contains;

template <Id I, Id H, Id... T> struct Contains<I, index_sequence<H, T...>> {
  static constexpr bool value =
      I == H || Contains<I, index_sequence<T...>>::value;
};

template <Id I> struct Contains<I, index_sequence<>> {
  static constexpr bool value = false;
};

} // namespace fn

template <class IT, class NT, int VO, typename Tries, typename Nodes>
struct FreeJoin;

// base case: no nodes
template <class IT, class NT, int VO, typename... Tries>
struct FreeJoin<IT, NT, VO, tuple<Tries...>, NodeList<>> {
  template <class Output>
  static void run(const tuple<Tries...> &tries, array<IT, VO> &tuple,
                  Value<NT> &value, Output &&output) {
    output(tuple, value);
  }
};

template <class IT, class NT, int VO, typename... Tries, typename... Nodes>
struct FreeJoin<IT, NT, VO, tuple<Tries...>, NodeList<Nodes...>> {
  static constexpr size_t num_tries = sizeof...(Tries);

  using AllTries = tuple<Tries...>;
  using AllTriesSeq = std::make_index_sequence<num_tries>;
  using Node = typename NodeList<Nodes...>::Head;
  using NodeSeq = std::make_index_sequence<Node::size()>;

  template <class T> using rpt = std::remove_pointer_t<T>;

  template <size_t Out, size_t In, size_t... Perm, size_t... Iota>
  static void project(array<IT, Out> &out, const array<IT, In> &in,
                      index_sequence<Perm...>, index_sequence<Iota...>) {
    static_assert(sizeof...(Perm) == sizeof...(Iota),
                  "permutation size mismatch");
    static_assert(sizeof...(Perm) <= Out, "output size too small");
    static_assert(sizeof...(Perm) <= In, "input size too small");
    // print
    // (..., (std::cout << Perm << " = in[" << Iota << "] = " << in[Iota]));
    // std::cout << std::endl;
    (..., (out[Perm] = in[Iota]));
  }

  template <size_t Out, size_t In, size_t... Perm>
  static void scatter(array<IT, Out> &out, const array<IT, In> &in,
                      index_sequence<Perm...> perm) {
    project(out, in, perm, std::make_index_sequence<perm.size()>());
  }

  template <size_t Out, size_t In, size_t... Perm>
  static void gather(array<IT, Out> &out, const array<IT, In> &in,
                     index_sequence<Perm...> perm) {
    project(out, in, std::make_index_sequence<perm.size()>(), perm);
  }

  // helper to merge the subtries into new_tries
  // new_tries[Node[0]] = subtries[0] ...
  template <typename OutSeq, typename InSeq> struct MergeTries;

  template <size_t... Out, size_t... In>
  struct MergeTries<index_sequence<Out...>, index_sequence<In...>> {
    template <typename AllTries, typename Subtries>
    static void run(AllTries &new_tries, Subtries &subtries) {
      (..., (get<Out>(new_tries) = get<In>(subtries)));
    }
  };

  template <class Output>
  static void run(const AllTries &tries, array<IT, VO> &tuple, Value<NT> &value,
                  Output &&output) {
    auto node_tries = fn::Map<Node, AllTries>::run(
        tries, [](auto trie, auto i) { return trie; });

    // constexpr size_t bound_vars =
    //     rpt<std::tuple_element_t<0, decltype(node_tries)>>::num_vars;

    // swap the minimal cover to the front
    // FIXME: this doesn't work if the covers have different variable
    // permutations... let's just assume they don't for now
    // size_t min_size = std::numeric_limits<size_t>::max();
    // fn::Apply<NodeSeq, decltype(node_tries)>::run(
    //     node_tries, [&](auto trie, auto i) {
    //       if constexpr (rpt<decltype(trie)>::num_vars == bound_vars) {
    //         size_t size = trie->getEstimatedSize();
    //         if (size < min_size) {
    //           min_size = size;
    //           std::swap(get<0>(node_tries), get<i>(node_tries));
    //         }
    //       }
    //     });

    auto cover = get<0>(node_tries);
    Value<NT> v;

    // iterate over the cover
    const size_t iter_size = cover->getIterSize();
    for (size_t item = 0; item < iter_size; ++item) {
      auto new_tries =
          fn::Map<AllTriesSeq, AllTries>::run(tries, [](auto trie, auto i) {
            // if i is in node, return nullptr to child trie, otherwise self
            if constexpr (fn::Contains<i, Node>::value) {
              return static_cast<typename rpt<decltype(trie)>::ChildTrie *>(
                  nullptr);
            } else {
              return trie;
            }
          });
      auto subtries = fn::Map<NodeSeq, decltype(node_tries)>::run(
          node_tries, [](auto trie, auto i) {
            return static_cast<typename rpt<decltype(trie)>::ChildTrie *>(
                nullptr);
          });
      // probe the remaining tries in the node
      if (fn::ApplyUntil<NodeSeq, decltype(node_tries)>::run(
              node_tries, [&](auto trie, auto i) -> bool {
                if constexpr (i == 0) {
                  // iterate over the cover
                  typename rpt<decltype(trie)>::Tuple cover_tuple;
                  trie->iter(item, cover_tuple, v);
                  // write cover tuple into the output tuple
                  scatter(tuple, cover_tuple,
                          typename rpt<decltype(trie)>::OutPerm());
                  if (trie->last_level) {
                    value.some *= v.some;
                    get<i>(subtries) = nullptr;
                  } else {
                    get<i>(subtries) = trie->lookup(cover_tuple);
                  }
                  return false;
                } else {
                  // probe the other tries
                  typename rpt<decltype(trie)>::Tuple key;
                  gather(key, tuple, typename rpt<decltype(trie)>::OutPerm());
                  auto subtrie = trie->lookup(key);
                  if (!subtrie) {
                    return true;
                  }
                  get<i>(subtries) = subtrie;
                  return false;
                }
              })) {
        continue;
      }
      // merge
      MergeTries<Node, NodeSeq>::run(new_tries, subtries);
      // recurse
      FreeJoin<IT, NT, VO, decltype(new_tries),
               typename NodeList<Nodes...>::Tail>::run(new_tries, tuple, value,
                                                       output);
    }
  }
};

} // namespace jit

struct Atom {
  Variable relation;
  std::vector<Variable> vars;

  Atom() = default;

  Atom(Variable relation, std::vector<Variable> vars)
      : relation(relation), vars(vars) {}
};

using JoinNode = std::vector<Atom>;
using FreeJoinPlan = std::vector<JoinNode>;

namespace detail {

using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

// Trying to pre-allocate stuff. Trying to outsmart the allocator. This
// probably just makes things worse, but I'm in too deep now!
template <class IT, class NT> struct FJAux {
  struct FrameAux {
    unique_ptr<Dim[]> trie_perm;
    unique_ptr<LazyTrie<IT, NT> *[]> all_tries;
    unique_ptr<LazyTrie<IT, NT> *[]> subtries;
  };

  vector<unique_ptr<Dim[]>> arrs;
  vector<unique_ptr<LazyTrie<IT, NT>>> tries;
  unique_ptr<IT[]> tuple;
  unique_ptr<IT[]> key;
  unique_ptr<FrameAux[]> frames;

  void alloc(Dim vars, Dim nodes, Dim relations) {
    tuple = std::make_unique<IT[]>(vars);
    key = std::make_unique<IT[]>(vars);
    frames = std::make_unique<FrameAux[]>(nodes + 1);
    for (Dim i = 0; i < nodes + 1; ++i) {
      frames[i].trie_perm = std::make_unique<Dim[]>(relations);
      frames[i].all_tries = std::make_unique<LazyTrie<IT, NT> *[]>(relations);
      frames[i].subtries = std::make_unique<LazyTrie<IT, NT> *[]>(relations);
    }
    for (size_t i = 0; i < tries.size(); ++i) {
      frames[0].all_tries[i] = tries[i].get();
    }
  }
};

struct FJ {
  vector<Variable> head_vars;
  // TODO(@altanh): try using (padded) flat arrays instead
  // offsets into all_tries; tries[i][0] is the cover trie for the i-th node
  vector<vector<size_t>> tries;
  // proj[i][j][k] is the index of the k-th variable of the j-th atom of the
  // i-th node in the output tuple
  vector<vector<vector<size_t>>> proj;
  // number of covers in each node
  vector<size_t> num_covers;

  template <class IT, class NT>
  static FJ make(const vector<unique_ptr<LazyTrie<IT, NT>>> &all_tries,
                 const FreeJoinPlan &plan) {
    FJ fj;

    // Collect head variables from covers; this is our output schema.
    // For now it is flattened (i.e. columnar) but one could feasibly return
    // it as in factorized form as a trie.
    unordered_map<Variable, size_t> var_to_idx;
    for (const auto &node : plan) {
      const auto &atom = node[0];
      for (const auto &var : atom.vars) {
        fj.head_vars.push_back(var);
        var_to_idx[var] = fj.head_vars.size() - 1;
      }
    }

    // For each node, find indices of tries
    for (const auto &node : plan) {
      vector<size_t> node_tries;
      size_t num_covers = 0;
      for (const auto &atom : node) {
        for (size_t i = 0; i < all_tries.size(); i++) {
          if (all_tries[i]->getRelation().name == atom.relation) {
            node_tries.push_back(i);
            break;
          }
        }
        if (atom.vars.size() == node[0].vars.size()) {
          num_covers++;
        }
      }
      fj.num_covers.push_back(num_covers);
      fj.tries.push_back(std::move(node_tries));
    }

    // For each node, find projection indices
    for (const auto &node : plan) {
      vector<vector<size_t>> node_proj;
      for (const auto &atom : node) {
        vector<size_t> atom_proj;
        for (const auto &var : atom.vars) {
          atom_proj.push_back(var_to_idx[var]);
        }
        node_proj.push_back(std::move(atom_proj));
      }
      fj.proj.push_back(std::move(node_proj));
    }

    // fj.print(all_tries);

    return fj;
  }

  template <class IT, class NT>
  void print(const vector<unique_ptr<LazyTrie<IT, NT>>> &all_tries) {
    for (size_t i = 0; i < tries.size(); ++i) {
      std::cout << "Node " << i << std::endl;
      for (size_t j = 0; j < tries[i].size(); ++j) {
        std::cout << "  Atom " << j << std::endl;
        std::cout << "    Relation: " << tries[i][j] << "->"
                  << all_tries[tries[i][j]]->getRelation().name << std::endl;
        std::cout << "    Vars: ";
        for (const auto &var : proj[i][j]) {
          std::cout << var << "->" << head_vars[var] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
};

template <class IT, class NT, class Output>
void FreeJoin(const FJ &plan, Output output, IT *tuple, Value<NT> value,
              size_t node, FJAux<IT, NT> *aux) {
  if (node == plan.tries.size()) {
    output(tuple, value);
    return;
  }

  const auto &all_tries = aux->frames[node].all_tries;

  // auto node_tries = plan.tries[node];
  // auto proj = plan.proj[node];
  const auto &node_tries = plan.tries[node];
  const auto &proj = plan.proj[node];

  // TODO: is this worth it, or should I just permute the tries in place?
  auto &trie_perm = aux->frames[node].trie_perm;
  for (Dim i = 0; i < node_tries.size(); ++i) {
    trie_perm[i] = i;
  }

  // TODO: there is an interesting tradeoff here
  // TODO: allow repeated binding of variables to recover binary joins?
  // covers are all_tries[node_tries[0:plan.num_covers[node]]
  {
    size_t cover_size = std::numeric_limits<size_t>::max();
    size_t cover_idx = 0;
    for (size_t i = 0; i < plan.num_covers[node]; ++i) {
      size_t size = all_tries[node_tries[i]]->getEstimatedSize();
      if (size < cover_size) {
        cover_size = size;
        cover_idx = i;
      }
    }
    std::swap(trie_perm[0], trie_perm[cover_idx]);
    // std::swap(node_tries[0], node_tries[cover_idx]);
    // std::swap(proj[0], proj[cover_idx]);

    // sort the rest by size
    // std::sort(trie_perm.begin() + 1, trie_perm.end(),
    //           [&all_tries](Dim a, Dim b) {
    //             return all_tries[a]->getMapSize() <
    //             all_tries[b]->getMapSize();
    //           });
  }

  // TODO: allocations
  // vector<LazyTrie<IT, NT> *> subtries(node_tries.size(), nullptr);
  // vector<LazyTrie<IT, NT> *> new_tries(all_tries);

  auto &subtries = aux->frames[node].subtries;
  auto &new_tries = aux->frames[node + 1].all_tries;
  for (size_t i = 0; i < aux->tries.size(); ++i) {
    new_tries[i] = all_tries[i];
  }

  // LazyTrie<IT, NT> *cover = all_tries[node_tries[0]];
  LazyTrie<IT, NT> *cover = all_tries[node_tries[trie_perm[0]]];
  // IT *cover_tuple = tuple + proj[0][0];
  IT *cover_tuple = tuple + proj[trie_perm[0]][0];
  Value<NT> v;

  const size_t cover_size = cover->getIterSize();
  for (size_t i = 0; i < cover_size; ++i) {
    cover->iter(i, cover_tuple, &v);
    // now tuple is filled upto cover
    if (cover->isLastCover()) {
      value.some *= v.some;
      subtries[0] = nullptr;
    } else {
      subtries[0] = cover->lookup(cover_tuple);
    }
    bool miss = false;
    for (size_t j = 1; j < node_tries.size(); ++j) {
      // LazyTrie<IT, NT> *trie = all_tries[node_tries[j]];
      LazyTrie<IT, NT> *trie = all_tries[node_tries[trie_perm[j]]];
      // project tuple down to keys for trie
      // const auto &trie_proj = proj[j];
      const auto &trie_proj = proj[trie_perm[j]];
      for (size_t k = 0; k < trie_proj.size(); ++k) {
        aux->key[k] = tuple[trie_proj[k]];
      }
      LazyTrie<IT, NT> *subtrie = trie->lookup(aux->key.get());
      if (!subtrie) {
        miss = true;
        break;
      }
      subtries[j] = subtrie;
    }
    if (miss) {
      continue;
    }
    // update all_tries
    for (size_t j = 0; j < node_tries.size(); ++j) {
      // new_tries[node_tries[j]] = subtries[j];
      new_tries[node_tries[trie_perm[j]]] = subtries[j];
    }
    FreeJoin(plan, output, tuple, value, node + 1, aux);
  }
}

template <class IT, class NT>
bool isValidPlan(const vector<Relation<IT, NT>> &relations,
                 const FreeJoinPlan &plan) {
  unordered_map<Variable, unordered_set<Variable>> rel_vars;
  for (const auto &rel : relations) {
    for (const auto &var : rel.vars) {
      rel_vars[rel.name].insert(var);
    }
  }
  for (const auto &node : plan) {
    unordered_set<Variable> seen_rels;
    unordered_set<Variable> cover_vars(node[0].vars.begin(),
                                       node[0].vars.end());
    for (const auto &atom : node) {
      // check that the same relation does not appear in multiple atoms
      if (seen_rels.count(atom.relation)) {
        std::cerr << "Relation " << atom.relation
                  << " appears multiple times in the same node" << std::endl;
        return false;
      }
      // check that all variables are bound at most once, and atom 0 is a
      // cover
      seen_rels.insert(atom.relation);
      for (const auto &var : atom.vars) {
        if (!cover_vars.count(var)) {
          std::cerr << "Variable " << var << " is not covered for "
                    << atom.relation << std::endl;
          return false;
        }
        if (!rel_vars[atom.relation].count(var)) {
          std::cerr << "Variable " << var << " is bound multiple times for "
                    << atom.relation << std::endl;
          return false;
        }
        rel_vars[atom.relation].erase(var);
      }
    }
  }
  // check that all variables are bound at least once
  for (const auto &rel : relations) {
    if (!rel_vars[rel.name].empty()) {
      std::cerr << "Relation " << rel.name << " has unbound variables:";
      for (const auto &var : rel_vars[rel.name]) {
        std::cerr << " " << var;
      }
      std::cerr << std::endl;
      return false;
    }
  }
  return true;
}

template <class IT, class NT>
void buildTries(vector<Relation<IT, NT>> &relations, const FreeJoinPlan &plan,
                FJAux<IT, NT> *aux) {
  // Compute trie schemas and cover information
  unordered_map<Variable, Schema> schemas;
  unordered_map<Variable, bool> is_cover;
  for (const auto &node : plan) {
    for (const auto &atom : node) {
      schemas[atom.relation].push_back(atom.vars);
      // check if atom is a cover; the final value will be true iff the last
      // subatom of a relation is also a cover of its node.
      is_cover[atom.relation] = atom.vars.size() == node[0].vars.size();
    }
  }
  // Create tries
  auto &tries = aux->tries;
  for (auto &rel : relations) {
    aux->arrs.emplace_back(nullptr);
    aux->arrs.emplace_back(nullptr);
    auto &nvars_out = aux->arrs[aux->arrs.size() - 2];
    auto &perm_out = aux->arrs[aux->arrs.size() - 1];
    tries.push_back(std::make_unique<LazyTrie<IT, NT>>(
        LazyTrie<IT, NT>::fromSchema(&rel, schemas[rel.name],
                                     is_cover[rel.name], &nvars_out,
                                     &perm_out)));
  }
}

} // namespace detail

template <class IT, class NT, class Output>
void FreeJoin(std::vector<Relation<IT, NT>> &relations,
              const FreeJoinPlan &plan, Output output) {
  // Validate plan
  if (!detail::isValidPlan(relations, plan)) {
    throw std::runtime_error("Invalid plan");
  }

  // Build tries
  detail::FJAux<IT, NT> aux;
  detail::buildTries(relations, plan, &aux);

  // Preprocess
  detail::FJ fj = detail::FJ::make(aux.tries, plan);

  // Allocate auxiliary memory
  aux.alloc(fj.head_vars.size(), plan.size(), relations.size());

  // Run the join
  detail::FreeJoin(fj, output, aux.tuple.get(), Value<NT>(), 0, &aux);
}

} // namespace space
