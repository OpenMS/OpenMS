#include <OpenDIALibrary.h>
#include <OpenDIAWorkflow.h>
#include <malloc.h>
#include <chrono>
#include <iostream>
using namespace OpenMS;
static double inuse() { auto m = mallinfo2(); return (m.uordblks + m.hblkhd) / 1048576.0; }
int main(int, char** argv) {
  double base = inuse();
  { auto l = OpenDIA::readAnyLibrary(argv[1]); std::cout << "search Library MB=" << inuse()-base << " tr=" << l.transitions.size() << " comp=" << l.compounds.size() << "\n";
    double b2 = inuse(); { auto t = l.traml_ids; std::cout << "traml_ids MB=" << inuse()-b2 << "\n"; }
    double b3 = inuse(); { auto c = l.compounds; std::cout << "compounds MB=" << inuse()-b3 << "\n"; } }
  { double b = inuse(); auto l = OpenDIA::readAnyLibrary(argv[2]); std::cout << "scoring Library MB=" << inuse()-b << "\n"; }
  double b = inuse(); auto a = OpenDIA::AnalysisLibrary::read(argv[1], argv[2]); double total = inuse()-b;
  double b4 = inuse(); { auto p = a.precursors; std::cout << "precursors MB=" << inuse()-b4 << "\n"; }
  double b5 = inuse(); { auto i = a.index; std::cout << "index MB=" << inuse()-b5 << "\n"; }
  double b6 = inuse(); { auto t = a.traml_ids; std::cout << "analysis traml_ids MB=" << inuse()-b6 << "\n"; }
  std::cout << "AnalysisLibrary total MB=" << total << " sizeof Precursor=" << sizeof(OpenDIA::AnalysisLibrary::Precursor) << "\n";
}
