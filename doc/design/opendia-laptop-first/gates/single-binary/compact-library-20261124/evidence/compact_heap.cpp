#include <OpenDIALibrary.h>
#include <OpenDIAWorkflow.h>
#include <malloc.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <unistd.h>
using namespace OpenMS;
static double heap() { auto m = mallinfo2(); return (m.uordblks + m.hblkhd) / 1048576.0; }
static double rss() { std::ifstream f("/proc/self/statm"); size_t total=0, resident=0; f >> total >> resident; return resident * double(sysconf(_SC_PAGESIZE)) / 1048576.0; }
template<class T> static void component(const char* name, const T& value) { double b=heap(); { auto copy=value; std::cout << name << " copy_heap_MiB=" << heap()-b << " size=" << copy.size() << "\n"; } }
int main(int argc, char** argv) {
  if (argc != 3) return 2;
  std::cout.setf(std::ios::unitbuf);
  double base=heap();
  { auto l=OpenDIA::readAnyLibrary(argv[1]); std::cout << "search Library heap_MiB=" << heap()-base << " tr=" << l.transitions.size() << " comp=" << l.compounds.size() << "\n"; component("search traml_ids", l.traml_ids); component("search compounds", l.compounds); component("search transitions", l.transitions); }
  { double b=heap(); auto l=OpenDIA::readAnyLibrary(argv[2]); std::cout << "scoring Library heap_MiB=" << heap()-b << " tr=" << l.transitions.size() << " comp=" << l.compounds.size() << "\n"; }
  double b=heap(); auto t=std::chrono::steady_clock::now(); auto a=OpenDIA::AnalysisLibrary::read(argv[1],argv[2]); double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count();
  std::cout << "AnalysisLibrary total heap_MiB=" << heap()-b << " retained_rss_MiB=" << rss() << " seconds=" << seconds << " sizeof_Precursor=" << sizeof(OpenDIA::AnalysisLibrary::Precursor) << " sizeof_Fragment=" << sizeof(OpenDIA::AnalysisLibrary::Fragment) << "\n";
  component("precursors", a.precursors); component("index", a.index); component("fragments", a.fragments); component("fragment_offset", a.fragment_offset); component("protein_sets", a.protein_sets); component("calibration_pairs", a.calibration_pairs); component("sibling_group_sizes", a.sibling_group_sizes);
}
