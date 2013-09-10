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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

///////////////////////////

START_TEST(OnDiscMSExperiment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

OnDiscMSExperiment<>* ptr = 0;
OnDiscMSExperiment<>* nullPointer = 0;
START_SECTION((OnDiscMSExperiment()))
{
	ptr = new OnDiscMSExperiment<>(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION((~OnDiscMSExperiment()))
{
	delete ptr;
}
END_SECTION

START_SECTION((OnDiscMSExperiment(const OnDiscMSExperiment& source)))
{
	OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> tmp2(tmp);
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getName(), tmp.getExperimentalSettings()->getInstrument().getName() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getVendor(), tmp.getExperimentalSettings()->getInstrument().getVendor() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getModel(), tmp.getExperimentalSettings()->getInstrument().getModel() )
  TEST_EQUAL(tmp2.getExperimentalSettings()->getInstrument().getMassAnalyzers().size(), tmp.getExperimentalSettings()->getInstrument().getMassAnalyzers().size() )
  TEST_EQUAL(tmp2.size(),tmp.size());
}
END_SECTION

START_SECTION((bool operator== (const OnDiscMSExperiment& rhs) const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> same(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));

  TEST_EQUAL(tmp==same, true);
  TEST_EQUAL((*tmp.getExperimentalSettings())==(*same.getExperimentalSettings()), true);
  TEST_EQUAL(tmp==failed, false);
}
END_SECTION

START_SECTION((bool operator!= (const OnDiscMSExperiment& rhs) const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> same(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));

  TEST_EQUAL(tmp!=same, false);
  TEST_EQUAL(tmp!=failed, true);
}
END_SECTION

START_SECTION((bool isSortedByRT() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.isSortedByRT(), true);
}
END_SECTION

START_SECTION((inline Size size() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.size(), 2);
  TEST_EQUAL(failed.size(), 0);
}
END_SECTION

START_SECTION((inline bool empty() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  TEST_EQUAL(failed.empty(), true);
}
END_SECTION

START_SECTION((inline Size getNrSpectra() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrSpectra(), 2);
  TEST_EQUAL(failed.getNrSpectra(), 0);
}
END_SECTION

START_SECTION((inline Size getNrChromatograms() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  OnDiscMSExperiment<> failed(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrChromatograms(), 1);
  TEST_EQUAL(failed.getNrChromatograms(), 0);
}
END_SECTION

START_SECTION((boost::shared_ptr<const ExperimentalSettings> getExperimentalSettings() const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  boost::shared_ptr<const ExperimentalSettings> settings = tmp.getExperimentalSettings();

  TEST_EQUAL(settings->getInstrument().getName(), "LTQ FT")
  TEST_EQUAL(settings->getInstrument().getMassAnalyzers().size(), 1)
}
END_SECTION

START_SECTION((inline MSSpectrum<PeakT>& operator[] (Size n) const))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  MSSpectrum<> s = tmp[0];
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);
}
END_SECTION

START_SECTION((MSSpectrum<PeakT> getSpectrum(Size id)))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.empty(), false);
  MSSpectrum<> s = tmp.getSpectrum(0);
  TEST_EQUAL(s.empty(), false);
  TEST_EQUAL(s.size(), 19914);
}
END_SECTION

START_SECTION((MSSpectrum<PeakT> getChromatogram(Size id)))
{
  OnDiscMSExperiment<> tmp(OPENMS_GET_TEST_DATA_PATH("IndexedmzMLFile_1.mzML"));
  TEST_EQUAL(tmp.getNrChromatograms(), 1);
  TEST_EQUAL(tmp.empty(), false);
  MSChromatogram<> c = tmp.getChromatogram(0);
  TEST_EQUAL(c.empty(), false);
  TEST_EQUAL(c.size(), 48);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

