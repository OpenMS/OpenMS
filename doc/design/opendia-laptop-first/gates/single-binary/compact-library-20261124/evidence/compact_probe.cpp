#include <OpenDIALibrary.h>
#include <OpenDIAWorkflow.h>
#include <malloc.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <unistd.h>
using namespace OpenMS;
static double heap() { auto m=mallinfo2(); return (m.uordblks+m.hblkhd)/1048576.0; }
static double rss() { std::ifstream f("/proc/self/statm"); size_t total=0,resident=0; f>>total>>resident; return resident*double(sysconf(_SC_PAGESIZE))/1048576.0; }
int main(int argc,char** argv) {
  if(argc<3) return 2;
  std::cout.setf(std::ios::unitbuf);
  double b=heap(); auto t=std::chrono::steady_clock::now();
  if(std::string(argv[1])=="analysis") {
    if(argc!=4) return 2;
    auto l=OpenDIA::AnalysisLibrary::read(argv[2],argv[3],[](const std::string& s){std::cout<<s<<"\n";});
    std::cout<<"analysis precursors="<<l.precursors.size()<<" fragments="<<l.fragments.size()<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<" heap_MiB="<<heap()-b<<" retained_rss_MiB="<<rss()<<"\n";
  } else {
    auto l=OpenDIA::readAnyLibrary(argv[2]);
    std::cout<<"read transitions="<<l.transitions.size()<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-t).count()<<" heap_MiB="<<heap()-b<<" retained_rss_MiB="<<rss()<<"\n";
  }
}
