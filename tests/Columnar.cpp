#include <space/format/Columnar.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>

using IT = uint64_t;
using NT = double;
using std::array;
using std::pair;
using std::vector;

using namespace space;

bool test_coo() {
  vector<array<IT, 3>> coo = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  auto data =
      Columnar<IT, NT>::fromCOO(coo, [](auto x) { return x[0] + x[1]; });

  if (data.getSize() != 3) {
    std::cerr << "size mismatch" << std::endl;
    std::cerr << "  expected: " << 3 << std::endl;
    std::cerr << "  got: " << data.getSize() << std::endl;
    return false;
  }

  if (data.getRank() != 3) {
    std::cerr << "rank mismatch" << std::endl;
    std::cerr << "  expected: " << 3 << std::endl;
    std::cerr << "  got: " << data.getRank() << std::endl;
    return false;
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      if (data.dim(j)[i] != coo[i][j]) {
        std::cerr << "dim mismatch" << std::endl;
        std::cerr << "  expected: " << coo[i][j] << std::endl;
        std::cerr << "  got: " << data.dim(j)[i] << std::endl;
        return false;
      }
    }
    if (*data.val(i) != coo[i][0] + coo[i][1]) {
      std::cerr << "val mismatch" << std::endl;
      std::cerr << "  expected: " << coo[i][0] + coo[i][1] << std::endl;
      std::cerr << "  got: " << data.val(i) << std::endl;
      return false;
    }
  }

  return true;
}

bool test_no_value() {
  vector<array<IT, 3>> coo = {
      {1, 2, 3},
      {4, 5, 6},
      {7, 8, 9},
  };
  auto data = Columnar<IT, void>::fromCOO(coo, [](auto x) {});

  if (data.getSize() != 3) {
    std::cerr << "size mismatch" << std::endl;
    std::cerr << "  expected: " << 3 << std::endl;
    std::cerr << "  got: " << data.getSize() << std::endl;
    return false;
  }

  if (data.getRank() != 3) {
    std::cerr << "rank mismatch" << std::endl;
    std::cerr << "  expected: " << 3 << std::endl;
    std::cerr << "  got: " << data.getRank() << std::endl;
    return false;
  }

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      if (data.dim(j)[i] != coo[i][j]) {
        std::cerr << "dim mismatch" << std::endl;
        std::cerr << "  expected: " << coo[i][j] << std::endl;
        std::cerr << "  got: " << data.dim(j)[i] << std::endl;
        return false;
      }
    }
  }

  return true;
}

int main(int argc, char **argv) {
  if (!test_coo()) {
    return 1;
  }
  if (!test_no_value()) {
    return 1;
  }
  return 0;
}
