#pragma once
#include <iostream>
#include <chrono>
#include <thread>

class Stopwatch {
private:
    std::chrono::steady_clock::time_point then{};
    std::chrono::milliseconds accumulator{};

public:
    void tic();
    void toc();
    void reset();
    std::chrono::milliseconds getElapsed() const;
};