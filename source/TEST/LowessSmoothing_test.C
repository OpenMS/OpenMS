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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LowessSmoothing* ptr = 0;
LowessSmoothing* null_ptr = 0;
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

std::vector<DoubleReal> x, y, out;

DoubleReal exp1[20] = {10, 20, 30, 40, 50, 60, 70.75057075, 81.91410776, 90.49938042, 93.67188213, 90.49938042, 81.91410776, 70.75057075, 60, 50, 40, 30, 20, 10, 7.638334409e-14};
DoubleReal exp2[20] = {4.940778184, 19.1953138, 32.45871201, 44.62566121, 55.59150285, 65.28588352, 73.78027456, 81.64413917, 87.38364167, 89.36964666, 87.38364167, 81.64413917, 72.86539444, 63.49165214, 53.94643243, 43.76172539, 32.89091229, 21.38760603, 9.323517923, -3.233540303};


for (DoubleReal i = 1.0; i <= 20.0; i += 1.0)
{
    x.push_back(i);
}

for (DoubleReal i = 1.0; i <= 10.0; i += 1.0)
{
    y.push_back(i*10);
}

for (DoubleReal i = 1.0; i <= 10.0; i += 1.0)
{
    y.push_back(100.0 - i*10);
}

y.push_back(10.0);

LowessSmoothing lowsmooth;
Param lowpar;
lowpar.setValue("window_size", 15);

START_SECTION(void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
    TEST_EXCEPTION(Exception::InvalidValue, lowsmooth.smoothData(x, y, out));

    y.pop_back();
    out.clear();

    lowsmooth.smoothData(x, y, out);

    for (Size i = 0; i < out.size(); ++i)
    {
        TEST_REAL_SIMILAR(out[i], exp1[i]);
    }

    out.clear();
    lowsmooth.setParameters(lowpar);
    lowsmooth.smoothData(x, y, out);

    for (Size i = 0; i < out.size(); ++i)
    {
        TEST_REAL_SIMILAR(out[i], exp2[i]);
    }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



