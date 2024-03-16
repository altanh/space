cmake_minimum_required(VERSION 3.10)


project(space)

# set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# don't omit frame pointers for better perfomance analysis
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

include_directories(include)
# print CXX flags for RelWithDebInfo
# message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")

# print active CXX flags
# message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

add_subdirectory(tests)
