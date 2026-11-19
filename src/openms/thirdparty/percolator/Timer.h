
#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>
#include <ctime>
#include <istream>
#include <string>

namespace OpenMS { namespace Internal { namespace Percolator {


/* Timer is for measuring wall clock time and CPU time of different code sections. Stop is used to make a checkpoint in time and reset for restarting the timer. */
class Timer
{

public:
    Timer();
    void reset();
    void stop();

    std::string getCPUTimeStr();
    std::string getWallTimeStr();
    std::string getStartTimeStr();

private:
    time_t wallEndTime, wallStartTime;
    clock_t cpuEndTime, cpuStartTime;
    std::string getTimeStr(double time, int numDecimals);
};

/* Memory use of the process, in KiB; -1 where it is not available. */
struct MemoryUsage {
    long long rssKb = -1;      // resident set size now (VmRSS)
    long long peakRssKb = -1;  // largest resident set size so far (VmHWM)
    long long heapKb = -1;     // allocated and not freed by malloc; can be
                               // above rssKb for untouched allocations
};

/* The resident set size and its peak from /proc/self/status (Linux), and
   the heap in use from mallinfo2 (glibc 2.33 or later). */
MemoryUsage getMemoryUsage();
/* Reads VmRSS and VmHWM from the text of a /proc/<pid>/status file. */
void parseProcStatus(std::istream& status, MemoryUsage& usage);
/* ", rss 988 MiB, peak 998 MiB, heap 950 MiB", leaving out what is not
   available. */
std::string formatMemoryUsage(const MemoryUsage& usage);

/* StageTimer reports the wall clock and CPU time spent in its scope, and the
   memory use at its end, to stderr when the verbosity is above 2. Used to see
   which processing stage dominates time and memory on large inputs. */
class StageTimer
{

public:
    explicit StageTimer(const char* name);
    ~StageTimer();

private:
    const char* name_;
    std::chrono::steady_clock::time_point wallStart_;
    clock_t cpuStart_;
};

}}}  // namespace OpenMS::Internal::Percolator

#endif

