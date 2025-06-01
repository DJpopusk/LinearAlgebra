#pragma once
#include <thread>

class ThreadPoolConfig {
public:
    static void setNumThreads(size_t num_threads);
    static size_t getNumThreads();
    
private:
    static size_t num_threads_;
};