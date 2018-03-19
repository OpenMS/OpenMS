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
#include <OpenMS/METADATA/IdentificationHit.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(IdentificationHit, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IdentificationHit* ptr = nullptr;
IdentificationHit* null_ptr = nullptr;
START_SECTION(IdentificationHit())
{
	ptr = new IdentificationHit();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IdentificationHit())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~IdentificationHit()))
{
  // TODO
}
END_SECTION

START_SECTION((IdentificationHit(const IdentificationHit &source)))
{
  // TODO
}
END_SECTION

START_SECTION((IdentificationHit& operator=(const IdentificationHit &source)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IdentificationHit &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IdentificationHit &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setId(const String &id)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getId() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCharge(Int charge)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCalculatedMassToCharge(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getCalculatedMassToCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setExperimentalMassToCharge(double mz)))
{
  // TODO
}
END_SECTION

START_SECTION((double getExperimentalMassToCharge() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPassThreshold(bool pass)))
{
  // TODO
}
END_SECTION

START_SECTION((bool getPassThreshold() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setRank(Int rank)))
{
  // TODO
}
END_SECTION

START_SECTION((Int getRank() const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



