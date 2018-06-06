// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

double targetFunction(double x)
{
  return 10 + 20*x + 40*x*x;
}

START_TEST(LowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LowessSmoothing* ptr = nullptr;
LowessSmoothing* null_ptr = nullptr;
START_SECTION(LowessSmoothing())
{
	ptr = new LowessSmoothing();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~LowessSmoothing())
{
	delete ptr;
}
END_SECTION

/////
std::vector<double> x, y, y_noisy, out;

//exact data
for (double i = 1.0; i <= 20.0; i += 1.0)
{
    x.push_back(i);
    y.push_back(targetFunction(i));
}

//noisy data
// make some noise

boost::random::mt19937 rnd_gen_;
for (Size i = 0; i < y.size(); ++i)
{
  boost::normal_distribution<float> udist (y.at(i), 0.05);
  y_noisy.push_back( udist(rnd_gen_) );
}

LowessSmoothing lowsmooth;
Param lowpar;
lowpar.setValue("window_size", 15);

TOLERANCE_RELATIVE(1.0004);
TOLERANCE_ABSOLUTE(0.07);
START_SECTION(void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
    y.push_back(-1.0);
    TEST_EXCEPTION(Exception::InvalidValue, lowsmooth.smoothData(x, y, out));
    y.pop_back();
    out.clear();

    lowsmooth.smoothData(x, y, out);
    Size idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

    out.clear();
    lowsmooth.setParameters(lowpar);
    lowsmooth.smoothData(x, y, out);
    idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

    out.clear();
    lowsmooth.smoothData(x, y_noisy, out);
    idx = 1;
    for (Size i = 0; i < out.size(); ++i, ++idx)
    {
        TEST_REAL_SIMILAR(out[i], targetFunction(idx));
    }

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



