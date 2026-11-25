

#include "Timer.h"

#include <fstream>
#include <iostream>
#include <sstream>

#if defined(__GLIBC__)
#include <malloc.h>
#endif

#include "Globals.h"

namespace OpenMS { namespace Internal { namespace Percolator {


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

std::string Timer::getStartTimeStr(){
    return ctime(&wallStartTime);
}

void parseProcStatus(std::istream& status, MemoryUsage& usage) {
    std::string line;
    while (std::getline(status, line)) {
        long long* field = NULL;
        if (line.compare(0, 6, "VmRSS:") == 0) {
            field = &usage.rssKb;
        } else if (line.compare(0, 6, "VmHWM:") == 0) {
            field = &usage.peakRssKb;
        } else {
            continue;
        }
        std::istringstream value(line.substr(6));
        long long kb = -1;
        if (value >> kb) {
            *field = kb;
        }
    }
}

MemoryUsage getMemoryUsage() {
    MemoryUsage usage;
#if defined(__linux__)
    std::ifstream status("/proc/self/status");
    if (status) {
        parseProcStatus(status, usage);
    }
#endif
#if defined(__GLIBC__) && \
    (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 33))
    // the chunks in use in all arenas plus the memory-mapped chunks
    struct mallinfo2 info = mallinfo2();
    usage.heapKb = static_cast<long long>((info.uordblks + info.hblkhd) / 1024u);
#endif
    return usage;
}

std::string formatMemoryUsage(const MemoryUsage& usage) {
    std::ostringstream out;
    if (usage.rssKb >= 0) {
        out << ", rss " << usage.rssKb / 1024 << " MiB";
    }
    if (usage.peakRssKb >= 0) {
        out << ", peak " << usage.peakRssKb / 1024 << " MiB";
    }
    if (usage.heapKb >= 0) {
        out << ", heap " << usage.heapKb / 1024 << " MiB";
    }
    return out.str();
}

StageTimer::StageTimer(const char* name)
    : name_(name), wallStart_(std::chrono::steady_clock::now()), cpuStart_(clock()) {}

StageTimer::~StageTimer()
{
    if (VERB > 2) {
        double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart_).count();
        double cpu = ((double)(clock() - cpuStart_)) / (double)CLOCKS_PER_SEC;
        std::cerr << "Stage " << name_ << ": " << wall << " s wall, " << cpu << " s cpu"
                  << formatMemoryUsage(getMemoryUsage()) << std::endl;
    }
}

}}}  // namespace OpenMS::Internal::Percolator
