#pragma once

#include <chrono>
struct Timer
{
    std::chrono::high_resolution_clock::time_point t1, t2;//varibles for time record
    std::chrono::duration<double> time_span;
    void start()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }
    void stop()
    {
        t2 = std::chrono::high_resolution_clock::now();
    }
    double GetRuntime()//return time in second
    {
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);//std::chrono::nanoseconds
        return time_span.count();
    }
};