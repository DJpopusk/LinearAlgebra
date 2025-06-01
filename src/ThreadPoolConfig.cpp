#include "ThreadPoolConfig.hpp"

size_t ThreadPoolConfig::num_threads_ = std::thread::hardware_concurrency();

void ThreadPoolConfig::setNumThreads(size_t num_threads) {
    num_threads_ = (num_threads > 0) ? num_threads : 1;
}

size_t ThreadPoolConfig::getNumThreads() {
    return num_threads_;
}