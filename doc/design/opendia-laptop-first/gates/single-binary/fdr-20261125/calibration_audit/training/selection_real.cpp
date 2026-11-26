// Real-table selection diagnostic. No biological or entrapment labels are read.
#define main synthetic_main
#include "selection_isolation.cpp"
#undef main
#include <fstream>
int main(int argc,char** argv){
  if(argc!=2)return 1;
  CrossFitScorer::Input in; std::vector<float> score;
  std::ifstream f(argv[1]);uint32_t precursor,pair,fold,decoy;float s;
  while(f>>precursor>>pair>>fold>>decoy>>s){in.precursor.push_back(precursor);in.pair.push_back(pair);in.fold.push_back(fold);in.decoy.push_back(decoy);score.push_back(s);}
  std::cout<<"rows="<<in.rows()<<'\n';
  for(unsigned held=0;held<3;++held){
    auto shifted=score;
    for(size_t r=0;r<in.rows();++r)if(in.fold[r]==held)shifted[r]+=(in.decoy[r]?0.5f:-0.5f);
    auto global=select(in,score,held,false), local=select(in,score,held,true);
    auto altered=select(in,shifted,held,false),isolated=select(in,shifted,held,true);
    std::cout<<"held="<<held<<" global_pos="<<global.pos.size()<<" local_pos="<<local.pos.size()<<" global_neg="<<global.neg.size()<<" local_neg="<<local.neg.size()<<" global_local_pos_diff="<<diff(global.pos,local.pos)<<" global_local_neg_diff="<<diff(global.neg,local.neg)<<" perturb_global_pos_diff="<<diff(global.pos,altered.pos)<<" perturb_global_neg_diff="<<diff(global.neg,altered.neg)<<" perturb_local_pos_diff="<<diff(local.pos,isolated.pos)<<" perturb_local_neg_diff="<<diff(local.neg,isolated.neg)<<'\n';
    if(diff(local.pos,isolated.pos)||diff(local.neg,isolated.neg))return 2;
  }
}
