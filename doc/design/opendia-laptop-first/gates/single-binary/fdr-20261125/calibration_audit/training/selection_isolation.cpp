// Diagnostic only. Uses the production pooled TDC and HistGBM without modifying OpenDIA.
#include <OpenDIAGBM.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
using namespace OpenMS::OpenDIA;
struct Selection { std::vector<uint32_t> pos, neg; };
Selection select(const CrossFitScorer::Input& in, const std::vector<float>& score, unsigned held, bool exclude_first)
{
  std::vector<uint32_t> rows;
  for (uint32_t r=0; r<in.rows(); ++r) if (!exclude_first || in.fold[r]!=held) rows.push_back(r);
  auto c=CrossFitScorer::pooledPairedTDC(in, score, rows);
  Selection s;
  for(size_t i=0;i<c.ranked.size();++i) if(c.q[i]<=0.01 && !in.decoy[c.ranked[i]] && in.fold[c.ranked[i]]!=held) s.pos.push_back(c.ranked[i]);
  std::vector<uint32_t> dw;
  for(auto r:c.winners) if(in.decoy[r]) dw.push_back(r);
  std::sort(dw.begin(),dw.end(),[&](auto a,auto b){return score[a]!=score[b]?score[a]>score[b]:a<b;});
  size_t take=std::min(dw.size(),std::max<size_t>(1,std::ceil(0.02*dw.size())));
  for(size_t i=0;i<take;++i) if(in.fold[dw[i]]!=held) s.neg.push_back(dw[i]);
  std::sort(s.pos.begin(),s.pos.end()); std::sort(s.neg.begin(),s.neg.end()); return s;
}
size_t diff(const std::vector<uint32_t>& a,const std::vector<uint32_t>& b)
{std::vector<uint32_t> d; std::set_symmetric_difference(a.begin(),a.end(),b.begin(),b.end(),std::back_inserter(d)); return d.size();}
std::vector<double> fit(const CrossFitScorer::Input& in,const Selection& s)
{
  std::vector<uint8_t> labels(in.rows()); std::vector<double> weights(in.rows(),1);
  std::vector<uint32_t> tr,va;
  for(auto r:s.pos){labels[r]=1;weights[r]=double(s.neg.size())/s.pos.size(); ((in.pair_hash[r]%100<10)?va:tr).push_back(r);}
  for(auto r:s.neg) ((in.pair_hash[r]%100<10)?va:tr).push_back(r);
  std::sort(tr.begin(),tr.end());std::sort(va.begin(),va.end());
  HistGBM model; HistGBM::Options o; o.max_iter=40; o.threads=2;
  if (!s.pos.empty() && !s.neg.empty() && !tr.empty()) model.fit(in.features.data(),in.dimension,tr,va,labels,weights,o);
  else model.init=0.0;
  std::vector<double> out;
  for(size_t r=0;r<in.rows();++r) if(in.fold[r]==0) out.push_back(model.predict(in.features.data()+r*in.dimension));
  return out;
}
int main()
{
  CrossFitScorer::Input in; in.dimension=4;std::vector<float> s;
  std::mt19937_64 rng(51);std::normal_distribution<float> normal(0,1);
  for(uint32_t p=0;p<15000;++p) for(unsigned d=0;d<2;++d){
    in.precursor.push_back(2*p+d);in.pair.push_back(p);in.decoy.push_back(d);in.fold.push_back(p%3);in.pair_hash.push_back(CrossFitScorer::fnv1a(std::to_string(p)));
    float x=normal(rng)+(!d && p%5?4.0f:0.0f);s.push_back(x);
    for(unsigned j=0;j<4;++j)in.features.push_back(j==0?x:normal(rng));
  }
  auto changed=s;
  // Only scores of fold 0 change. Every other fold's scores/features/labels remain fixed.
  for(size_t r=0;r<in.rows();++r) if(in.fold[r]==0) changed[r]+=(in.decoy[r]?0.5f:-0.5f);
  bool witnessed=false;
  for(bool isolated:{false,true}){
    auto a=select(in,s,0,isolated),b=select(in,changed,0,isolated);
    auto pd=diff(a.pos,b.pos),nd=diff(a.neg,b.neg);
    auto pa=fit(in,a),pb=fit(in,b);size_t changed_probes=0;double max_delta=0;
    for(size_t i=0;i<pa.size();++i){double d=std::abs(pa[i]-pb[i]);max_delta=std::max(max_delta,d);changed_probes+=d>1e-12;}
    std::cout<<(isolated?"exclude_before_selection":"production_global_selection")<<" positives="<<a.pos.size()<<"/"<<b.pos.size()<<" negatives="<<a.neg.size()<<"/"<<b.neg.size()<<" changed_positive_rows="<<pd<<" changed_negative_rows="<<nd<<" changed_fixed_probe_predictions="<<changed_probes<<" max_delta="<<max_delta<<'\n';
    if(!isolated)witnessed=pd+nd>0 && changed_probes>0;
    else if(pd||nd||changed_probes)return 2;
  }
  return witnessed?0:3;
}
