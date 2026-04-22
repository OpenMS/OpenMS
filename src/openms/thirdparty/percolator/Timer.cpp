

#include "Timer.h"

#include <iostream>
#include <sstream>

Timer::Timer(){
    reset();
}

void Timer::reset()
{
    cpuStartTime = cpuEndTime = clock();
    time(&wallStartTime);
    time(&wallEndTime);
}
void Timer::stop()
{
    cpuEndTime = clock();
    time(&wallEndTime);
}

std::string Timer::getCPUTimeStr()
{
    double cpuTime = ((double)(cpuEndTime - cpuStartTime)) / (double)CLOCKS_PER_SEC;
    return getTimeStr(cpuTime, 4);
}

std::string Timer::getWallTimeStr()
{
    return getTimeStr(difftime(wallEndTime, wallStartTime), 0);
}

std::string Timer::getTimeStr(double time, int numDecimals)
{
    std::stringstream timeStr;
    timeStr.precision(numDecimals);
    timeStr << std::fixed << time;
    return timeStr.str();
}

char* Timer::getStartTimeStr(){
    return ctime(&wallStartTime);
}

