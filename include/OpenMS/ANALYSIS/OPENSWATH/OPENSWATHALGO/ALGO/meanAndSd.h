/*
 * meanAndSd.h
 *
 *  Created on: Jul 24, 2012
 *      Author: witold
 */

#ifndef MEANANDSD_H_
#define MEANANDSD_H_
#include <cmath>

namespace OpenSwath
{
	class mean_and_stddev
	{
			double m_, q_;
			unsigned long c_;
		public:
			typedef double argument_type, result_type;
			mean_and_stddev()
			: m_(0.0), q_(0.0), c_(0u)
			{
			}

			void operator ()(double sample)
			{
				double const delta = sample - m_;
				m_ += delta / ++c_;
				q_ += delta * (sample - m_);
			}

			double sample_variance() const
			{
				return (c_ > 1u) ? (q_ / (c_ - 1)) : 0;
			}
			double standard_variance() const
			{
				return (c_>1u) ? (q_ / c_) : 0;
			}
			double sample_stddev() const
			{
				return std::sqrt(sample_variance());
			}
			double standard_stddev() const
			{
				return std::sqrt(standard_variance());
			}

			double mean() const
			{
				return m_;
			}
			unsigned long count() const
			{
				return c_;
			}
			double variance() const
			{
				return sample_variance();
			}
			double stddev() const
			{
				return sample_stddev();
			}
			double operator ()() const
			{
				return stddev();
			}
	};



}

#endif /* MEANANDSD_H_ */
