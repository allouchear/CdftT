#ifndef CDFTT_TIMER_H_INCLUDED
#define CDFTT_TIMER_H_INCLUDED

#include <sys/time.h>

class Timer{
    timeval timerStart;
    timeval timerStop;
    public:
    void init()
    {
        gettimeofday(&timerStart, NULL);
    }
    Timer(){
        init();
    }
    double get()
    {
        gettimeofday(&timerStop, NULL);
        timeval timerElapsed;
        timersub(&timerStop, &timerStart, &timerElapsed);
        double t=timerElapsed.tv_sec*1000.0+timerElapsed.tv_usec/1000.0;
        timerStart=timerStop;
        return t;
    }
    
};

#endif //CDFTT_TIMER_H_INCLUDED