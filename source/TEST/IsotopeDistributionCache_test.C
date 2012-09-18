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
// $Maintainer: Lars Nilse $
// $Authors: Bastian Blank $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/FILTERING/DATAREDUCTION/IsotopeDistributionCache.h>

using namespace OpenMS;

START_TEST(IsotopeDistributionCache, "$Id$")

START_SECTION(IsotopeDistributionCache(DoubleReal max_mass, DoubleReal mass_window_width, DoubleReal intensity_percentage=0, DoubleReal intensity_percentage_optional=0))
  IsotopeDistributionCache c(100, 1);
END_SECTION 

START_SECTION(const TheoreticalIsotopePattern& getIsotopeDistribution(DoubleReal mass) const)
  IsotopeDistributionCache c(1000, 10);
  const IsotopeDistributionCache::TheoreticalIsotopePattern &p(c.getIsotopeDistribution(500));
  TEST_REAL_SIMILAR(p.intensity[0], 1);
  TEST_REAL_SIMILAR(p.intensity[1], 0.2668);
  TEST_REAL_SIMILAR(p.intensity[2], 0.04864);
  TEST_REAL_SIMILAR(p.intensity[3], 0.006652);
  TEST_EQUAL(&p == &c.getIsotopeDistribution(509.9), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(510.0), true);
  TEST_EQUAL(&p != &c.getIsotopeDistribution(499.9), true);
END_SECTION

END_TEST

