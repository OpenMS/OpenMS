#include <OpenDIALibrary.h>
#include <OpenDIAWorkflow.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <chrono>
#include <iostream>
using namespace OpenMS;
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static std::string mem() { size_t a=0,b=0; SysInfo::getProcessMemoryConsumption(a); SysInfo::getProcessPeakMemoryConsumption(b); return " rss_mb="+std::to_string(a/1024)+" peak_mb="+std::to_string(b/1024); }
int main(int argc, char** argv) {
  std::string mode = argv[1];
  double t = now();
  if (mode == "read") { auto l = OpenDIA::readAnyLibrary(argv[2], [](const std::string& s){std::cout<<s<<"\n";}); std::cout << "read " << l.transitions.size() << " s=" << now()-t << mem() << "\n"; }
  else if (mode == "analysis") { auto l = OpenDIA::AnalysisLibrary::read(argv[2], argv[3], [](const std::string& s){std::cout<<s<<"\n";}); std::cout << "analysis " << l.precursors.size() << " s=" << now()-t << mem() << "\n";
    t = now(); OpenDIA::Library copy = l.search; std::cout << "copy s=" << now()-t << mem() << "\n"; }
}
