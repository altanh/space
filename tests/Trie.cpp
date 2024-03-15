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
  const size_t size = trie.getSize();
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

int main(int argc, char **argv) {
  auto coo = vector<array<IT, 3>>{
      {1, 2, 3},
      {1, 2, 6},
      {7, 8, 9},
  };
  auto data = std::make_shared<Columnar<IT, NT>>(
      Columnar<IT, NT>::fromCOO(coo, [](auto x) { return x[0] + x[1]; }));
  Relation<IT, NT> A = {"A", {"x", "y", "z"}, data};

  Schema s1;
  s1.push_back({"x", "y"});
  s1.push_back({"z"});
  LazyTrie<IT, NT> t1(A, s1);
  printSchema(t1.getSchema());
  printTrie(t1, 0, t1.getSchema().size() - 1);
  std::cout << std::endl;

  Schema s2;
  s2.push_back({"x"});
  s2.push_back({"y", "z"});
  s2.push_back({});
  LazyTrie<IT, NT> t2(A, s2);
  printSchema(t2.getSchema());
  printTrie(t2, 0, t2.getSchema().size() - 1);
  std::cout << std::endl;

  Schema s3;
  s3.push_back({"x"});
  s3.push_back({"y"});
  s3.push_back({"z"});
  LazyTrie<IT, NT> t3(A, s3);
  printSchema(t3.getSchema());
  printTrie(t3, 0, t3.getSchema().size() - 1);
  std::cout << std::endl;

  Schema s4;
  s4.push_back({"x", "y", "z"});
  LazyTrie<IT, NT> t4(A, s4);
  printSchema(t4.getSchema());
  printTrie(t4, 0, t4.getSchema().size() - 1);
  std::cout << std::endl;

  return 0;
}