// Species-free mixed-null stress test, not a proteomics validation or a new decoy design.
// Null target and decoy features are iid. Known true targets have a shift in feature 0.
#include <OpenDIAGBM.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
using namespace OpenMS::OpenDIA;
int main(int argc,char**argv){
  unsigned seed=argc>1?std::stoul(argv[1]):0;
  CrossFitScorer::Input in;in.dimension=6;std::vector<float> score;std::vector<bool> truth;
  std::mt19937_64 rng(seed+910);std::normal_distribution<float> normal(0,1);
  for(uint32_t p=0;p<12000;++p)for(unsigned d=0;d<2;++d){
    in.precursor.push_back(2*p+d);in.pair.push_back(p);in.decoy.push_back(d);in.fold.push_back(p%3);in.pair_hash.push_back(CrossFitScorer::fnv1a(std::to_string(p)));
    bool signal=!d && p%5!=0;truth.push_back(signal);
    float x=normal(rng)+(signal?3.0f:0);score.push_back(x);
    for(unsigned j=0;j<6;++j)in.features.push_back(j==0?x:normal(rng));
  }
  std::vector<uint32_t> rows;for(uint32_t r=0;r<in.rows();++r)rows.push_back(r);
  auto report=[&](const char*tag,const auto&c){size_t n=0,e=0,dec=0;
    for(size_t i=0;i<c.ranked.size();++i)if(c.q[i]<=0.01){auto r=c.ranked[i];if(in.decoy[r])++dec;else{++n;e+=!truth[r];}}
    std::cout<<seed<<'\t'<<tag<<'\t'<<n<<'\t'<<e<<'\t'<<dec<<'\t'<<(n?double(e)/n:0)<<'\n';};
  report("seed",CrossFitScorer::pooledPairedTDC(in,score,rows));
  auto result=gbmRescore(in,score,0.02,seed,2);
  report("gbm",CrossFitScorer::pooledPairedTDC(in,result.score,rows));
}
