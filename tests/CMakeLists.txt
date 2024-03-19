# helper function to add a test
function(add_test name)
    add_executable(test-${name} ${name}.cpp ${HEADERS})
    target_compile_options(test-${name} PRIVATE -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -fno-omit-frame-pointer)
    target_link_libraries(test-${name} PRIVATE "${GRAPHBLAS_LIB}")
endfunction()

add_test(HashTable)
add_test(Columnar)
add_test(Trie)
add_test(FreeJoin)
