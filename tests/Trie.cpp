#include <space/format/Trie.h>

#include <array>
#include <cstdint>
#include <memory>
#include <tuple>
#include <vector>

using namespace space;
using IT = uint64_t;
using NT = double;
using std::array;
using std::tuple;
using std::vector;

void printSchema(const Schema &schema) {
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
}

void printTrie(LazyTrie<IT, NT> &trie, size_t level, size_t max_level) {
  const size_t size = trie.getIterSize();
  std::unique_ptr<IT[]> key = std::make_unique<IT[]>(trie.getNumVars());
  for (size_t i = 0; i < size; ++i) {
    Value<NT> value;
    trie.iter<IT, NT>(i, key.get(), &value);
    for (size_t j = 0; j < level; ++j) {
      std::cout << "    ";
    }
    std::cout << "(";
    if (trie.getNumVars() > 0) {
      std::cout << key[0];
    }
    for (size_t j = 1; j < trie.getNumVars(); ++j) {
      std::cout << ", " << key[j];
    }
    if (level == max_level) {
      std::cout << ") -> " << value.some << std::endl;
    } else {
      std::cout << ")" << std::endl;
      LazyTrie<IT, NT> *child = trie.lookup(key.get());
      if (!child) {
        std::cerr << "key not found: " << key[0] << ", " << key[1] << std::endl;
        return;
      }
      printTrie(*child, level + 1, max_level);
    }
  }
}

int main() {
  auto coo = vector<array<IT, 3>>{
      {1, 2, 3},
      {1, 2, 6},
      {7, 8, 9},
  };
  auto data =
      Columnar<IT, NT>::fromCOO(coo, [](auto x) { return x[0] + x[1]; });
  Relation<IT, NT> A = {"A", {"x", "y", "z"}, &data};

  std::unique_ptr<Dim[]> nvars;
  std::unique_ptr<Dim[]> perm;

  Schema s1;
  s1.push_back({"x", "y"});
  s1.push_back({"z"});
  auto t1 = LazyTrie<IT, NT>::fromSchema(&A, s1, true, &nvars, &perm);
  t1.printSchema();
  printTrie(t1, 0, t1.getNumLevels() - 1);
  std::cout << std::endl;

  Schema s2;
  s2.push_back({"x"});
  s2.push_back({"y", "z"});
  auto t2 = LazyTrie<IT, NT>::fromSchema(&A, s2, true, &nvars, &perm);
  t2.printSchema();
  printTrie(t2, 0, t2.getNumLevels() - 1);
  std::cout << std::endl;

  Schema s3;
  s3.push_back({"x"});
  s3.push_back({"y"});
  s3.push_back({"z"});
  auto t3 = LazyTrie<IT, NT>::fromSchema(&A, s3, true, &nvars, &perm);
  t3.printSchema();
  printTrie(t3, 0, t3.getNumLevels() - 1);
  std::cout << std::endl;

  Schema s4;
  s4.push_back({"x", "y", "z"});
  auto t4 = LazyTrie<IT, NT>::fromSchema(&A, s4, true, &nvars, &perm);
  t4.printSchema();
  printTrie(t4, 0, t4.getNumLevels() - 1);
  std::cout << std::endl;

  // [x, y], [z]
  StaticLazyTrie<IT, NT, true, true, std::index_sequence<2, 1>,
                 std::index_sequence<0, 1, 2>>
      st(&data);
  st.print();

  // [z], [y], [x]
  StaticLazyTrie<IT, NT, true, true, std::index_sequence<1, 1, 1>,
                 std::index_sequence<2, 1, 0>>
      st2(&data);
  st2.print();

  return 0;
}