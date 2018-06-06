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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIdentification* ptr = nullptr;
CompNovoIdentification* nullPointer = nullptr;
START_SECTION(CompNovoIdentification())
{
  ptr = new CompNovoIdentification();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~CompNovoIdentification())
{
  delete ptr;
}
END_SECTION

START_SECTION((CompNovoIdentification(const CompNovoIdentification& source)))
  CompNovoIdentification cni;
  Param p(cni.getParameters());
  p.setValue("fragment_mass_tolerance", 0.5);
  cni.setParameters(p);
  TEST_EQUAL(CompNovoIdentification(cni).getParameters() == p, true)
END_SECTION


START_SECTION((void getIdentifications(std::vector< PeptideIdentification > &ids, const PeakMap &exp)))
{
   TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_losses", "true");
  tsg_param.setValue("add_isotopes", "true");
  tsg.setParameters(tsg_param);

  PeakSpectrum rspec;
  tsg.getSpectrum(rspec, AASequence::fromString("DFPIANGER"), 1, 1);

  PeakSpectrum spec;
  for (Size i = 0; i != rspec.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec[i].getMZ());
    p.setIntensity(rspec[i].getIntensity());
    spec.push_back(p);
  }

  PeakSpectrum rspec_ETD;

  tsg_param.setValue("add_b_ions", "false");
  tsg_param.setValue("add_y_ions", "false");
  tsg_param.setValue("add_z_ions", "true");
  tsg.setParameters(tsg_param);
  tsg.getSpectrum(rspec_ETD, AASequence::fromString("DFPIANGER"), 1, 1);

  tsg_param.setValue("add_z_ions", "false");
  tsg_param.setValue("add_precursor_peaks", "true");
  tsg.setParameters(tsg_param);
  tsg.getSpectrum(rspec_ETD, AASequence::fromString("DFPIANGER"), 2, 2);

  PeakSpectrum spec_ETD;
  for (Size i = 0; i != rspec_ETD.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec_ETD[i].getMZ());
    p.setIntensity(rspec_ETD[i].getIntensity());
    spec_ETD.push_back(p);
  }

  Precursor prec;
  prec.setMZ((AASequence::fromString("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);
  spec_ETD.setPrecursors(precs);

  PeakMap exp;
  exp.addSpectrum(spec);
  exp.addSpectrum(spec_ETD);

  vector<PeptideIdentification> ids;
  CompNovoIdentification cni;
  Param cni_param(cni.getParameters());
  cni.setParameters(cni_param);
  cni.getIdentifications(ids, exp);
  TEST_EQUAL(ids.size(), 1)
  TEST_EQUAL(ids.begin()->getHits().size() > 0, true)
  // After mass correction for b1 ions (#1440) a different peptide scored best.
  TEST_EQUAL(ids.begin()->getHits().begin()->getSequence() == AASequence::fromString("DFPDALGQR"), true)
}
END_SECTION

START_SECTION((void getIdentification(PeptideIdentification &id, const PeakSpectrum &CID_spec, const PeakSpectrum &ETD_spec)))
{
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_losses", "true");
  tsg_param.setValue("add_isotopes", "true");
  tsg.setParameters(tsg_param);

  PeakSpectrum rspec;
  tsg.getSpectrum(rspec, AASequence::fromString("DFPIANGER"), 1, 1);

  PeakSpectrum spec;
  for (Size i = 0; i != rspec.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec[i].getMZ());
    p.setIntensity(rspec[i].getIntensity());
    spec.push_back(p);
  }

  PeakSpectrum rspec_ETD;

  tsg_param.setValue("add_b_ions", "false");
  tsg_param.setValue("add_y_ions", "false");
  tsg_param.setValue("add_z_ions", "true");
  tsg.setParameters(tsg_param);
  tsg.getSpectrum(rspec_ETD, AASequence::fromString("DFPIANGER"), 1, 1);

  tsg_param.setValue("add_z_ions", "false");
  tsg_param.setValue("add_precursor_peaks", "true");
  tsg.setParameters(tsg_param);
  tsg.getSpectrum(rspec_ETD, AASequence::fromString("DFPIANGER"), 2, 2);

  PeakSpectrum spec_ETD;
  for (Size i = 0; i != rspec_ETD.size(); ++i)
  {
    Peak1D p;
    p.setMZ(rspec_ETD[i].getMZ());
    p.setIntensity(rspec_ETD[i].getIntensity());
    spec_ETD.push_back(p);
  }

  Precursor prec;
  prec.setMZ((AASequence::fromString("DFPLANGER").getMonoWeight() + 2.0 * Constants::PROTON_MASS_U) / 2.0);
  prec.setCharge(2);
  vector<Precursor> precs;
  precs.push_back(prec);
  spec.setPrecursors(precs);
  spec_ETD.setPrecursors(precs);

  PeptideIdentification id;
  CompNovoIdentification cni;
  Param cni_param(cni.getParameters());
  cni.setParameters(cni_param);
  cni.getIdentification(id, spec, spec_ETD);
  TEST_EQUAL(id.getHits().size() > 0, true)
  // After mass correction for b1 ions (#1440) a different peptide scored best.
  std::cout << id.getHits().begin()->getSequence() << std::endl;
  TEST_EQUAL(id.getHits().begin()->getSequence() == AASequence::fromString("DFPDALGQR"), true)

}
END_SECTION

START_SECTION((CompNovoIdentification& operator=(const CompNovoIdentification &source)))
{
  CompNovoIdentification cni;
  Param p(cni.getParameters());
  p.setValue("fragment_mass_tolerance", 0.5);
  cni.setParameters(p);
  CompNovoIdentification cni2;
  cni2 = cni;
  TEST_EQUAL(cni2.getParameters() == p, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



