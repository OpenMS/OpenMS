#include <torch/torch.h>
#include <cstdio>
#include <vector>
namespace {
#include "helpers.inc"
#include "coelution.inc"
};
std::vector<float> input(int64_t n, int64_t c, int64_t l, uint64_t s)
{
  std::vector<float> v(n * c * l);
  for (auto& x : v) { s = s * 6364136223846793005ULL + 1442695040888963407ULL; x = static_cast<float>((s >> 33) % 1000); }
  return v;
}
void print(const char* name, const torch::Tensor& t)
{
  auto c = t.contiguous();
  for (int64_t i = 0; i < c.size(0); ++i)
  {
    std::printf("%s", name);
    for (int64_t j = 0; j < c.size(1); ++j) { std::printf("\t%.6f", c[i][j].item<float>()); }
    std::printf("\n");
  }
}
int main()
{
  const int64_t n = 3, l = 32;
  auto a = input(n, 6, l, 7), b = input(n, 3, l, 11);
  auto ms2 = torch::from_blob(a.data(), {n, 6, l}).clone(), ms1 = torch::from_blob(b.data(), {n, 3, l}).clone();
  Coelution c6(6), c3(3);
  // Zero the learned similarity embedding to match the Rust zero-weight VarBuilder.
  torch::NoGradGuard g;
  for (auto& p : c6->parameters()) { p.zero_(); }
  for (auto& p : c3->parameters()) { p.zero_(); }
  print("coe6", c6(ms2));
  print("coe3", c3(ms1));
  print("ms12", crossModal(ms1, ms2));
}
