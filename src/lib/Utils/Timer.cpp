//#include <iostream>

#include <Utils/Timer.hpp>


Timer::Timer()
{
    init();
}

void Timer::init()
{
    gettimeofday(&timerStart, nullptr);
}

double Timer::get()
{
    gettimeofday(&timerStop, nullptr);

    timeval timerElapsed;
    timersub(&timerStop, &timerStart, &timerElapsed);

    double t = timerElapsed.tv_sec * 1000.0 + timerElapsed.tv_usec / 1000.0;
    timerStart = timerStop;

    return t;
}
