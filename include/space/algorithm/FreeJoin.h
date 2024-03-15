#pragma once

#include <space/format/Trie.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace space {

struct Atom {
  Variable relation;
  std::vector<Variable> vars;

  Atom(Variable relation, std::vector<Variable> vars)
      : relation(relation), vars(vars) {}
};

using JoinNode = std::vector<Atom>;
using FreeJoinPlan = std::vector<JoinNode>;

namespace detail {

using std::unordered_map;
using std::vector;

struct FJ {
  vector<Variable> head_vars;
  // offsets into all_tries; tries[i][0] is the cover trie for the i-th node
  vector<vector<size_t>> tries;
  // proj[i][j][k] is the index of the k-th variable of the j-th atom of the
  // i-th node in the output tuple
  vector<vector<vector<size_t>>> proj;

  template <class IT, class NT>
  void print(const vector<LazyTrie<IT, NT> *> &all_tries) {
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

template <class IT, class NT>
FJ preprocessFJ(const std::vector<LazyTrie<IT, NT> *> &all_tries,
                const FreeJoinPlan &plan) {
  FJ fj;

  // Collect head variables from covers; this is our output schema.
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
    for (const auto &atom : node) {
      for (size_t i = 0; i < all_tries.size(); i++) {
        if (all_tries[i]->getRelation().name == atom.relation) {
          node_tries.push_back(i);
          break;
        }
      }
    }
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

  fj.print(all_tries);

  return fj;
}

template <class IT, class NT, class Output>
void FreeJoin(const vector<LazyTrie<IT, NT> *> &all_tries, const FJ &plan,
              Output output, IT *tuple, Value<NT> *value, size_t node,
              IT *workspace) {
  if (node == plan.tries.size()) {
    output(tuple, value);
    return;
  }

  const auto &node_tries = plan.tries[node];
  // TODO: choose "smallest" cover using vector estimate per paper
  // TODO: there is an interesting tradeoff here
  LazyTrie<IT, NT> *cover = all_tries[node_tries[0]];
  const size_t cover_size = cover->getSize();
  IT *cover_tuple = tuple + plan.proj[node][0][0];
  Value<NT> v;

  vector<LazyTrie<IT, NT> *> subtries(node_tries.size(), nullptr);
  vector<LazyTrie<IT, NT> *> new_tries(all_tries);

  for (size_t i = 0; i < cover_size; ++i) {
    cover->iter(i, cover_tuple, &v);
    // now tuple is filled upto cover
    if (cover->isSuffix()) {
      value->some *= v.some;
      subtries[0] = nullptr;
    } else {
      subtries[0] = cover->lookup(cover_tuple);
    }
    bool miss = false;
    for (size_t j = 1; j < node_tries.size(); ++j) {
      LazyTrie<IT, NT> *trie = all_tries[node_tries[j]];
      // project tuple down to keys for trie
      for (size_t k = 0; k < plan.proj[node][j].size(); ++k) {
        workspace[k] = tuple[plan.proj[node][j][k]];
      }
      LazyTrie<IT, NT> *subtrie = trie->lookup(workspace);
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
      new_tries[node_tries[j]] = subtries[j];
    }
    FreeJoin(new_tries, plan, output, tuple, value, node + 1, workspace);
  }
}

} // namespace detail

template <class IT, class NT, class Output>
void FreeJoin(std::vector<LazyTrie<IT, NT> *> tries, const FreeJoinPlan &plan,
              Output output) {
  // Preprocess
  detail::FJ fj = detail::preprocessFJ(tries, plan);
  // TODO: NT *tuple = new NT[batch_size * nvars];
  IT *tuple = new IT[fj.head_vars.size()];
  IT *workspace = new IT[fj.head_vars.size()];
  Value<NT> value;
  detail::FreeJoin(tries, fj, output, tuple, &value, 0, workspace);
  delete[] tuple;
  delete[] workspace;
}

template <class IT, class NT, class Output>
void FreeJoin(const std::vector<Relation<IT, NT>> &relations,
              const FreeJoinPlan &plan, Output output) {
  // Compute trie schemas
  std::unordered_map<Variable, Schema> schemas;
  for (const auto &node : plan) {
    for (const auto &atom : node) {
      schemas[atom.relation].push_back(atom.vars);
    }
  }
  // Create tries
  std::vector<LazyTrie<IT, NT> *> tries;
  for (const auto &rel : relations) {
    std::unordered_map<Variable, size_t> var_to_idx;
    for (size_t i = 0; i < rel.vars.size(); ++i) {
      var_to_idx[rel.vars[i]] = i;
    }
    std::vector<Variable> transposed_vars;
    std::unique_ptr<size_t[]> perm =
        std::make_unique<size_t[]>(rel.vars.size());
    size_t i = 0;
    for (const auto &vars : schemas[rel.name]) {
      for (const auto &var : vars) {
        transposed_vars.push_back(var);
        perm[i++] = var_to_idx[var];
      }
    }
    // print permutation
    for (size_t i = 0; i < rel.vars.size(); ++i) {
      std::cout << perm[i] << " ";
    }
    std::cout << std::endl;
    // use transposed relation
    Relation<IT, NT> transposed = {
        rel.name, transposed_vars,
        std::make_shared<Columnar<IT, NT>>(*rel.data, std::move(perm))};
    tries.push_back(new LazyTrie<IT, NT>(transposed, schemas[rel.name]));
  }

  // Preprocess
  detail::FJ fj = detail::preprocessFJ(tries, plan);
  // TODO: NT *tuple = new NT[batch_size * nvars];
  IT *tuple = new IT[fj.head_vars.size()];
  IT *workspace = new IT[fj.head_vars.size()];
  Value<NT> value = {1.0};
  detail::FreeJoin(tries, fj, output, tuple, &value, 0, workspace);
  delete[] tuple;
  delete[] workspace;

  for (auto trie : tries) {
    delete trie;
  }
}

} // namespace space
