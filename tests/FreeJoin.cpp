#include <array>
#include <space/algorithm/FreeJoin.h>
#include <vector>

using namespace space;
using std::array;
using std::vector;
using IT = uint64_t;
using NT = double;

int main(int argc, char **argv) {
  auto coo = vector<array<IT, 2>>{{1, 2}, {2, 3}, {3, 1}};
  auto data = std::make_shared<Columnar<IT, NT>>(
      Columnar<IT, NT>::fromCOO(coo, [](auto x) { return 1.0; }));

  auto R = Relation<IT, NT>{"R", {"x", "y"}, data};
  auto S = Relation<IT, NT>{"S", {"y", "z"}, data};
  auto T = Relation<IT, NT>{"T", {"z", "x"}, data};

  FreeJoinPlan plan;
  plan.push_back({Atom("R", {"x"}), Atom("T", {"x"})});
  plan.push_back({Atom("R", {"y"}), Atom("S", {"y"})});
  plan.push_back({Atom("S", {"z"}), Atom("T", {"z"})});

  std::vector<Relation<IT, NT>> all_relations{R, S, T};
  FreeJoin(all_relations, plan, [](IT *idx, Value<NT> *value) {
    for (size_t i = 0; i < 3; i++) {
      std::cout << idx[i] << " ";
    }
    std::cout << " -> " << value->some << std::endl;
  });

  // auto R_schema = Schema{{"x"}, {"y"}};
  // auto S_schema = Schema{{"y"}, {"z"}};
  // auto T_schema = Schema{{"x"}, {"z"}};

  // auto R_trie = LazyTrie<IT, NT>(R, R_schema);
  // auto S_trie = LazyTrie<IT, NT>(S, S_schema);
  // auto T_trie = LazyTrie<IT, NT>(T, T_schema);

  // auto all_tries = vector<LazyTrie<IT, NT> *>{&R_trie, &S_trie, &T_trie};

  // FreeJoin(all_tries, plan, [](IT *idx, Value<NT> *value) {
  //   for (size_t i = 0; i < 3; i++) {
  //     std::cout << idx[i] << " ";
  //   }
  //   std::cout << " -> " << value->some << std::endl;
  // });

  return 0;
}