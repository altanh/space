#include <space/format/Trie.h>

#include <array>
#include <cstdint>
#include <tuple>
#include <vector>

using namespace space;
using IT = uint64_t;
using NT = double;
using std::array;
using std::tuple;
using std::vector;

int main(int argc, char **argv) {
  auto coo = vector<array<IT, 3>>{
      {1, 2, 3},
      {1, 2, 6},
      {7, 8, 9},
  };
  auto data = std::make_shared<Columnar<IT, NT>>(
      Columnar<IT, NT>::fromCOO(coo, [](auto x) { return x[0] + x[1]; }));

  Relation<IT, NT> A = {"A", {"x", "y", "z"}, data};
  Schema schema;
  schema.push_back({"x", "y"});
  schema.push_back({"z"});

  LazyTrie<IT, NT> trie(A, schema, true);

  const size_t size = trie.getSize();
  std::cout << "size: " << size << std::endl;

  for (size_t i = 0; i < size; ++i) {
    IT key[2];
    double value = 1.0;
    trie.iter<NT>(i, key, &value);
    std::cout << "L0: " << key[0] << ", " << key[1] << std::endl;
    LazyTrie<IT, NT> *child = trie.lookup(key);
    if (!child) {
      std::cerr << "key not found: " << key[0] << ", " << key[1] << std::endl;
      return 1;
    }
    const size_t child_size = child->getSize();
    for (size_t j = 0; j < child_size; ++j) {
      IT child_key[1];
      child->iter<NT>(j, child_key, &value);
      std::cout << "\tL1: " << child_key[0] << " -> " << value << std::endl;
    }
  }

  return 0;
}