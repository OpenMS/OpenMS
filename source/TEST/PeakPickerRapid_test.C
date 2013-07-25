// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerRapid.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerRapid, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerRapid* ptr = 0;
PeakPickerRapid* null_ptr = 0;
START_SECTION(PeakPickerRapid())
{
    ptr = new PeakPickerRapid();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~PeakPickerRapid()))
{
  delete ptr;
}
END_SECTION

Peak1D p1, p2, p3;
p1.setMZ(100.5);
p1.setIntensity(0.3520653);

p2.setMZ(101.0);
p2.setIntensity(0.3989423);

p3.setMZ(101.6);
p3.setIntensity(0.3332246);

PeakPickerRapid ppr;

START_SECTION((template < typename PeakType > bool computeTPG(const PeakType &p1, const PeakType &p2, const PeakType &p3, DoubleReal &mu, DoubleReal &sigma, DoubleReal &area, DoubleReal &height) const ))
{
    DoubleReal m, s, A, h;

    ppr.computeTPG(p1, p2, p3, m, s, A, h);

    std::cout << "results: " << m << " " << s << " " << A << " " << h << std::endl;

}
END_SECTION

START_SECTION((template < typename PeakType > void pick(const MSSpectrum< PeakType > &cinput, MSSpectrum< PeakType > &output)))
{
  // TODO
}
END_SECTION

START_SECTION((template < typename PeakType > void pickExperiment(MSExperiment< PeakType > &input, MSExperiment< PeakType > &output)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



