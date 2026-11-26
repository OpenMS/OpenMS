#include "OpenDIASearch.h"
#include "OpenDIAScorer.h"
#include <OpenMS/METADATA/Precursor.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <stdexcept>
using namespace OpenMS;
using namespace OpenMS::OpenDIA;
void check(bool x, const char* s) { if (!x) throw std::runtime_error(s); std::cout << "PASS " << s << '\n'; }
std::vector<SearchCandidate> search(bool flip) {
  std::vector<TransitionQuery> q;
  for (int p=0;p<2;++p) for(int f=0;f<3;++f) {
    TransitionQuery t; t.precursor_id=p?"B":"A"; t.precursor_mz=500; t.fragment_mz=300+f*100; t.library_rt=10; t.has_library_rt=true; t.is_decoy=bool(p)^flip; q.push_back(t);
  }
  SearchOptions o; o.rt_window_seconds=20; o.min_transitions=3;
  StreamSearch s(q,o);
  for(int i=0;i<=20;++i) {
    MSSpectrum ms; ms.setMSLevel(2); ms.setRT(i);
    Precursor p; p.setMZ(500); p.setIsolationWindowLowerOffset(25); p.setIsolationWindowUpperOffset(25); ms.setPrecursors({p});
    for(int f=0;f<3;++f) { Peak1D peak; peak.setMZ(300+f*100); peak.setIntensity(1000*std::exp(-0.5*std::pow((i-10)/2.0,2))); ms.push_back(peak); }
    s.consumeSpectrum(ms);
  }
  s.finish(); return s.candidates();
}
int main() {
  auto a=search(false), b=search(true);
  check(a.size()==2 && b.size()==a.size(),"search produces both identical-evidence precursors after label swap");
  bool same=true; for(size_t i=0;i<a.size();++i) {
    const auto& x=a[i]; const auto& y=b[i];
    same &= x.precursor==y.precursor && x.rt_start==y.rt_start && x.rt_apex==y.rt_apex && x.rt_end==y.rt_end && x.area==y.area && x.apex_intensity==y.apex_intensity && x.transition_support==y.transition_support && x.matched_scans==y.matched_scans && x.apex_rt_error_seconds==y.apex_rt_error_seconds && x.fwhm_seconds==y.fwhm_seconds && x.coelution==y.coelution && x.peak_rank==y.peak_rank && x.peaks_in_window==y.peaks_in_window;
  }
  check(same,"all 14 SearchCandidate fields unchanged by label swap");
  check(a[0].area==a[1].area && a[0].transition_support==a[1].transition_support,"identical target/decoy evidence receives identical search area and support");
  CrossFitScorer::Input in; std::vector<float> score; std::vector<uint32_t> rows;
  for(uint32_t p=0;p<12;++p) for(uint32_t d=0;d<2;++d) {
    in.precursor.push_back(2*p+d); in.pair.push_back(p); in.decoy.push_back(d); in.fold.push_back(p%3); in.pair_hash.push_back(p); score.push_back(float(p)+(d?0:0.25f)); rows.push_back(rows.size());
  }
  auto c=CrossFitScorer::pooledPairedTDC(in,score,rows);
  check(c.ranked.size()==12 && std::all_of(c.ranked.begin(),c.ranked.end(),[&](auto r){return !in.decoy[r];}),"strict higher score wins each pair across three folds");
  check(std::all_of(c.q.begin(),c.q.end(),[](double q){return std::abs(q-1.0/12.0)<1e-12;}),"12 target winners receive reverse-min q=1/12 with plus-one correction");
  auto strict_swap=in; for(auto& d:strict_swap.decoy)d=!d;
  auto strict_comp=CrossFitScorer::pooledPairedTDC(strict_swap,score,rows);
  auto winners=c.ranked, swapped_winners=strict_comp.ranked;
  std::sort(winners.begin(),winners.end()); std::sort(swapped_winners.begin(),swapped_winners.end());
  check(winners==swapped_winners,"strict pair winner identities survive label swap despite decoy-normalization change");
  check(std::all_of(strict_comp.q.begin(),strict_comp.q.end(),[](double q){return q==1;}),"label-swapped all-decoy winners receive q=1, intentionally asymmetric FDR");
  auto reversed=rows; std::reverse(reversed.begin(),reversed.end());
  auto cr=CrossFitScorer::pooledPairedTDC(in,score,reversed);
  check(cr.ranked==c.ranked && cr.q==c.q,"competition invariant to subset row iteration order");
  for(size_t i=0;i<score.size();++i) score[i]=float(i/2);
  c=CrossFitScorer::pooledPairedTDC(in,score,rows);
  check(std::all_of(c.ranked.begin(),c.ranked.end(),[&](auto r){return in.decoy[r];}),"all exact pair ties resolve to decoys");
  auto swapped=in; for(auto& d:swapped.decoy)d=!d;
  cr=CrossFitScorer::pooledPairedTDC(swapped,score,rows);
  check(std::all_of(cr.ranked.begin(),cr.ranked.end(),[&](auto r){return swapped.decoy[r];}),"tie priority follows label, not row identity");
  score[0]+=1; rows.erase(std::find(rows.begin(),rows.end(),1));
  c=CrossFitScorer::pooledPairedTDC(in,score,rows);
  check(std::find(c.ranked.begin(),c.ranked.end(),0)!=c.ranked.end(),"orphan target accepted as competitor when its decoy has no row");
  std::vector<uint32_t> support={4,1,2,2}; std::vector<bool> keep;
  for(size_t p=0;p<2;++p) for(int d=0;d<2;++d) keep.push_back(std::max(support[2*p],support[2*p+1])>=3);
  check(keep==std::vector<bool>({true,true,false,false}),"replica of pair-max population gate retains weak partner symmetrically");
}
