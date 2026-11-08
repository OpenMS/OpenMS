#include <torch/torch.h>
#include <cstdio>
namespace {
#include "helpers.inc"
#include "branch.inc"
#include "coelution.inc"
}
int64_t count(const torch::nn::Module& m) { int64_t n = 0; for (const auto& p : m.parameters()) n += p.numel(); return n; }
int main() { for (int64_t c : {6, 3}) { Branch b(c, 32); Coelution h(c); std::printf("branch_c%ld\t%ld\tcoelution_c%ld\t%ld\n", c, count(*b), c, count(*h)); } }
