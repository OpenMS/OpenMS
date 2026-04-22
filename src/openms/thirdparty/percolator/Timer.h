
#ifndef TIMER_H_
#define TIMER_H_

#include <string>
#include <ctime>

/* Timer is for measuring wall clock time and CPU time of different code sections. Stop is used to make a checkpoint in time and reset for restarting the timer. */
class Timer
{

public:
    Timer();
    void reset();
    void stop();

    std::string getCPUTimeStr();
    std::string getWallTimeStr();
    char* getStartTimeStr();

private:
    time_t wallEndTime, wallStartTime;
    clock_t cpuEndTime, cpuStartTime;
    std::string getTimeStr(double time, int numDecimals);
};

#endif

