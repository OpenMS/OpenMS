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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentificationCID, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIdentificationCID* ptr = 0;
CompNovoIdentificationCID* nullPointer = 0;
START_SECTION(CompNovoIdentificationCID())
{
	ptr = new CompNovoIdentificationCID();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~CompNovoIdentificationCID())
{
	delete ptr;
}
END_SECTION

START_SECTION((CompNovoIdentificationCID(const CompNovoIdentificationCID& source)))
	CompNovoIdentificationCID cnis;
	Param p(cnis.getParameters());
	p.setValue("precursor_mass_tolerance", 3.0);
	cnis.setParameters(p);
	TEST_EQUAL(CompNovoIdentificationCID(cnis).getParameters() == p, true)
END_SECTION

START_SECTION((CompNovoIdentificationCID& operator = (const CompNovoIdentificationCID& source)))
  CompNovoIdentificationCID cnis;
  Param p(cnis.getParameters());
  p.setValue("precursor_mass_tolerance", 3.0);
  cnis.setParameters(p);
	CompNovoIdentificationCID cnis2;
	cnis2 = cnis;
  TEST_EQUAL(cnis2.getParameters() == p, true)
END_SECTION

START_SECTION((void getIdentifications(std::vector<PeptideIdentification>& ids, const PeakMap& exp)))
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_losses", "true");
  tsg_param.setValue("add_isotopes", "true");
  tsg.setParameters(tsg_param);

  RichPeakSpectrum rspec;
  tsg.getSpectrum(rspec, AASequence("DFPIANGER"));

  PeakSpectrum spec;
  for (Size i = 0; i != rspec.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec[i].getMZ());
    p.setIntensity(rspec[i].getIntensity());
    spec.push_back(p);
  }

  Precursor prec;
  prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);

  vector<PeptideIdentification> ids;
  CompNovoIdentificationCID cni;
  Param cni_param(cni.getParameters());
  cni_param.setValue("precursor_mass_tolerance", 0.3);
  cni.setParameters(cni_param);
	PeakMap exp;
	exp.push_back(spec);
  cni.getIdentifications(ids, exp);
  TEST_EQUAL(ids.size(), 1)
  TEST_EQUAL(ids.begin()->getHits().size() > 0, true)
  TEST_EQUAL(ids.begin()->getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)
END_SECTION

START_SECTION((void getIdentification(PeptideIdentification& id, const PeakSpectrum& CID_spec)))
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_losses", "true");
  tsg_param.setValue("add_isotopes", "true");
  tsg.setParameters(tsg_param);

  RichPeakSpectrum rspec;
  tsg.getSpectrum(rspec, AASequence("DFPIANGER"));

  PeakSpectrum spec;
  for (Size i = 0; i != rspec.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec[i].getMZ());
    p.setIntensity(rspec[i].getIntensity());
    spec.push_back(p);
  }

	Precursor prec;
	prec.setMZ((AASequence("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
	prec.setCharge(2);
	vector<Precursor> precs;
	precs.push_back(prec);
	spec.setPrecursors(precs);	

	PeptideIdentification id;
	CompNovoIdentificationCID cni;
	Param cni_param(cni.getParameters());
	cni_param.setValue("precursor_mass_tolerance", 0.3);
	cni.setParameters(cni_param);
	cni.getIdentification(id, spec);
	TEST_EQUAL(id.getHits().size() > 0, true)
	TEST_EQUAL(id.getHits().begin()->getSequence() == AASequence("DFPLANGER"), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



