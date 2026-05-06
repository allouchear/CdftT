#ifndef CDFTT_TIMER_HPP_INCLUDED
#define CDFTT_TIMER_HPP_INCLUDED

#include <sys/time.h>

class Timer
{
    timeval timerStart;
    timeval timerStop;

    public:
        Timer();

        void init();

        double get();
};

#endif //CDFTT_TIMER_HPP_INCLUDED