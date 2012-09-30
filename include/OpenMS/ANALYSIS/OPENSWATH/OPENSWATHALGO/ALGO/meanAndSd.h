// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_MEANANDSD_H_
#define OPENMS_ANALYSIS_OPENSWATH_OPENSWATHALGO_MEANANDSD_H_
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

#endif 
