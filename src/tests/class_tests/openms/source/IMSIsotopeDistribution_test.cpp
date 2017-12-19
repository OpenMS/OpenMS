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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSIsotopeDistribution.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSIsotopeDistribution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSIsotopeDistribution* ptr = nullptr;
IMSIsotopeDistribution* null_ptr = nullptr;
START_SECTION(IMSIsotopeDistribution())
{
	ptr = new IMSIsotopeDistribution();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSIsotopeDistribution())
{
	delete ptr;
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(nominal_mass_type nominalMass=0)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(mass_type mass)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(const peaks_container &peaks, nominal_mass_type nominalMass=0)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((size_type size() const ))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator=(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const IMSIsotopeDistribution &distribution) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const IMSIsotopeDistribution &distribution) const ))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator*=(const IMSIsotopeDistribution &distribution)))
{
  // TODO
}
END_SECTION

START_SECTION((IMSIsotopeDistribution& operator*=(unsigned int pow)))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getMass(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((abundance_type getAbundance(size_type i) const ))
{
  // TODO
}
END_SECTION

START_SECTION((mass_type getAverageMass() const ))
{
  // TODO
}
END_SECTION

START_SECTION((nominal_mass_type getNominalMass() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setNominalMass(nominal_mass_type nominalMass)))
{
  // TODO
}
END_SECTION

START_SECTION((masses_container getMasses() const ))
{
  // TODO
}
END_SECTION

START_SECTION((abundances_container getAbundances() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void normalize()))
{
  // TODO
}
END_SECTION

START_SECTION((bool empty() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ims::IMSIsotopeDistribution::Peak] Peak(mass_type mass=0.0, abundance_type abundance=0.0)))
{
  // TODO
}
END_SECTION

START_SECTION(([ims::IMSIsotopeDistribution::Peak] bool operator==(const Peak &peak) const ))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



