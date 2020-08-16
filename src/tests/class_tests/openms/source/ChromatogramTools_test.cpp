// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/KERNEL/StandardTypes.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ChromatogramTools, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ChromatogramTools* ptr = nullptr;
ChromatogramTools* nullPointer = nullptr;
START_SECTION(ChromatogramTools())
{
	ptr = new ChromatogramTools();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~ChromatogramTools())
{
	delete ptr;
}
END_SECTION

START_SECTION((ChromatogramTools(const ChromatogramTools &)))
{
  ChromatogramTools tmp;
	ChromatogramTools tmp2(tmp);
	NOT_TESTABLE
}
END_SECTION

START_SECTION(template <typename ExperimentType> void convertChromatogramsToSpectra(ExperimentType& exp))
{
  PeakMap exp;
	MSChromatogram chrom1, chrom2;
	chrom1.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
	Precursor pre1, pre2;
	pre1.setMZ(100.1);
	pre2.setMZ(100.2);

	Product pro1, pro2;
	pro1.setMZ(200.1);
	pro2.setMZ(200.2);

	chrom1.setPrecursor(pre1);
	chrom1.setProduct(pro1);

	chrom2.setPrecursor(pre2);
	chrom2.setProduct(pro2);		

	chrom2.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
	ChromatogramPeak peak1, peak2, peak3;
	peak1.setRT(0.1);
	peak2.setRT(0.2);
	peak3.setRT(0.3);
	chrom1.push_back(peak1);
	chrom1.push_back(peak2);

	chrom2.push_back(peak2);
	chrom2.push_back(peak2);

	exp.addChromatogram(chrom1);
	exp.addChromatogram(chrom2);

	TEST_EQUAL(exp.size(), 0)
	TEST_EQUAL(exp.getChromatograms().size(), 2)
	ChromatogramTools().convertChromatogramsToSpectra(exp);
	TEST_EQUAL(exp.size(), 4)
	TEST_EQUAL(exp.getChromatograms().size(), 0)
	TEST_REAL_SIMILAR(exp[0][0].getMZ(), 200.1)

	TEST_EQUAL(exp[0].getPrecursors().size(), 1)
	TEST_REAL_SIMILAR(exp[0].getPrecursors().begin()->getMZ(), 100.1)


}
END_SECTION

START_SECTION(template <typename ExperimentType> void convertSpectraToChromatograms(ExperimentType& exp, bool remove_spectra = false))
{
  PeakSpectrum spec1, spec2, spec3, spec4, spec5;
	spec1.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
	spec2.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
	spec3.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
	spec4.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
	spec5.getInstrumentSettings().setScanMode(InstrumentSettings::MASSSPECTRUM);

	Precursor prec1, prec2;
	prec1.setMZ(500.1);
	prec2.setMZ(500.2);
	

	Peak1D p;
	p.setMZ(100.1);
	p.setIntensity(20000000);
	spec1.push_back(p);
	spec1.setRT(0.1);
	spec1.getPrecursors().push_back(prec1);

	p.setMZ(100.2);
	p.setIntensity(30000000);
	spec2.push_back(p);
	spec2.setRT(0.3);
	spec2.getPrecursors().push_back(prec2);

	p.setMZ(100.1);
	p.setIntensity(40000000);
	spec3.push_back(p);
	spec3.setRT(0.4);
	spec3.getPrecursors().push_back(prec1);

	p.setMZ(100.2);
	p.setIntensity(50000000);
	spec4.push_back(p);
	spec4.setRT(0.5);
	spec4.getPrecursors().push_back(prec2);

	PeakMap exp;
	exp.addSpectrum(spec1);
	exp.addSpectrum(spec2);
	exp.addSpectrum(spec3);
	exp.addSpectrum(spec4);
	exp.addSpectrum(spec5);

	PeakMap exp2 = exp;

	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(exp.getChromatograms().size(), 0)
	ChromatogramTools().convertSpectraToChromatograms(exp);
	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(exp.getChromatograms().size(), 2)

	TEST_EQUAL(exp2.size(), 5)
	TEST_EQUAL(exp2.getChromatograms().size(), 0)
	ChromatogramTools().convertSpectraToChromatograms(exp2, true);
	TEST_EQUAL(exp2.size(), 1)
	TEST_EQUAL(exp2.getChromatograms().size(), 2)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



