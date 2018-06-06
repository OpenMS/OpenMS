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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFileOptions, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFileOptions* ptr = nullptr;
FeatureFileOptions* null_ptr = nullptr;
START_SECTION(FeatureFileOptions())
{
	ptr = new FeatureFileOptions();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFileOptions())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setIntensityRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasIntensityRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getIntensityRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setLoadConvexHull(bool convex)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getLoadConvexHull() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setLoadSubordinates(bool sub)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getLoadSubordinates() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMetadataOnly(bool only)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getMetadataOnly() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setSizeOnly(bool only)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getSizeOnly() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRTRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasRTRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getRTRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setMZRange(const DRange< 1 > &range)))
{
  // TODO
}
END_SECTION

START_SECTION((bool hasMZRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((const DRange<1>& getMZRange() const ))
{
  // TODO
}
END_SECTION

START_SECTION((~FeatureFileOptions()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



