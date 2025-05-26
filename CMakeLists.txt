cmake_minimum_required(VERSION 3.10)
project(Linear_alg_Threads)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra)

include_directories(include)
set(SOURCES
        src/main.cpp
        src/Matrix.cpp
)

set(HEADERS
        include/Matrix.hpp
        include/Exceptions.hpp
)

add_executable(linear_solver1 ${SOURCES})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(linear_solver1 PRIVATE Threads::Threads)