// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MathFunctions.h>
#include <boost/random/mersenne_twister.hpp> // for mt19937_64
#include <boost/random/uniform_int.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/complement.hpp>

namespace OpenMS::Math
{
struct RandomShuffler::Impl
{ boost::mt19937_64 rng; };

RandomShuffler::RandomShuffler(): impl_(std::make_unique<Impl>())
{
}
RandomShuffler::RandomShuffler(int value): RandomShuffler()
{ impl_->rng.seed(value); }
RandomShuffler::RandomShuffler(const RandomShuffler& other): impl_(other.impl_ ? std::make_unique<Impl>(*other.impl_) : std::make_unique<Impl>())
{
}
RandomShuffler::RandomShuffler(RandomShuffler&& other) noexcept = default;
RandomShuffler& RandomShuffler::operator=(const RandomShuffler& other)
{
  if (this != &other) impl_ = other.impl_ ? std::make_unique<Impl>(*other.impl_) : std::make_unique<Impl>();
  return *this;
}
RandomShuffler& RandomShuffler::operator=(RandomShuffler&& other) noexcept = default;
RandomShuffler::~RandomShuffler() = default;
void RandomShuffler::seed(uint64_t value)
{
  if (! impl_) impl_ = std::make_unique<Impl>();
  impl_->rng.seed(value);
}
std::ptrdiff_t RandomShuffler::randomIndex_(std::ptrdiff_t upper)
{
  if (! impl_) impl_ = std::make_unique<Impl>();
  boost::uniform_int<std::ptrdiff_t> distribution(0, upper);
  return distribution(impl_->rng);
}

double log_binomial_coef(unsigned n, unsigned k)
{
  // Handle edge cases for improved numerical stability
  if (k > n) { throw std::invalid_argument("k cannot be greater than n in binomial coefficient"); }

  if (k == 0 || k == n)
  {
    return 0.0; // log(1) = 0
  }

  // Use symmetry to minimize computation for large k
  if (k > n / 2) { k = n - k; }

  return boost::math::lgamma(n + 1.0) - boost::math::lgamma(k + 1.0) - boost::math::lgamma(n - k + 1.0);
}

double binomial_cdf_complement(unsigned N, unsigned n, double p)
{
  if (p < 0.0 || p > 1.0) { throw std::invalid_argument("Probability p must be between 0 and 1"); }
  if (n > N) { throw std::invalid_argument("n cannot be greater than N"); }

  if (n == 0) return 1.0; // P(X ≥ 0) = 1
  if (p == 0.0) return (n == 0) ? 1.0 : 0.0;
  if (p == 1.0) return 1.0; // all mass at N

  const boost::math::binomial_distribution<double> dist(N, p);
  return boost::math::cdf(boost::math::complement(dist, n - 1));
}
} // namespace OpenMS::Math
