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
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIonScoring* ptr = nullptr;
CompNovoIonScoring* nullPointer = nullptr;
START_SECTION(CompNovoIonScoring())
{
	ptr = new CompNovoIonScoring();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((CompNovoIonScoring(const CompNovoIonScoring &source)))
{
  CompNovoIonScoring cnis;
	Param p(cnis.getParameters());
	p.setValue("fragment_mass_tolerance", 0.5);
	cnis.setParameters(p);
	TEST_EQUAL(CompNovoIonScoring(cnis).getParameters() == p, true)
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoring()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void scoreSpectra(Map< double, IonScore > &CID_ion_scores, PeakSpectrum &CID_spec, PeakSpectrum &ETD_spec, double precursor_weight, Size charge)))
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

	Map<double, CompNovoIonScoringBase::IonScore> ion_scores;
	CompNovoIonScoring cnis;
  cnis.scoreSpectra(ion_scores, spec, spec_ETD, 1018.48, 1);

  for (Map<double, CompNovoIonScoringBase::IonScore>::ConstIterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
  {
/*
y1 175.118952187571
y2 304.161545285171
y3 361.183009010571
y4 475.225936461371
y5 546.263050250571
y6 659.347114231171
y7 756.399878084171
y8 903.468292000971

b1 117.042044532471
b2 263.102633417371
b3 360.155397270371
b4 473.239461250971
b5 544.276575040171
b6 658.319502490971
b7 715.340966216371
b8 844.383559313971

*/
    cerr << it->first << " " << it->second.score << endl;
    if (fabs(it->first - 903.468292000971) < 0.001 ||
        fabs(it->first - 756.399878084171) < 0.001 ||
        fabs(it->first - 659.347114231171) < 0.001 ||
				fabs(it->first - 659.328) < 0.001 ||
        fabs(it->first - 546.263050250571) < 0.001 ||
        fabs(it->first - 475.225936461371) < 0.001 ||
        fabs(it->first - 361.183009010571) < 0.001 ||
				fabs(it->first - 361.164) < 0.001 ||
        fabs(it->first - 304.161545285171) < 0.001 ||
        fabs(it->first - 175.118952187571) < 0.001 ||
        fabs(it->first - 263.102633417371) < 0.001 ||
        fabs(it->first - 360.155397270371) < 0.001 ||
        fabs(it->first - 473.239461250971) < 0.001 ||
        fabs(it->first - 544.276575040171) < 0.001 ||
        fabs(it->first - 658.319502490971) < 0.001 ||
        fabs(it->first - 715.340966216371) < 0.001 ||
        fabs(it->first - 844.383559313971) < 0.001 ||
        //After introducing mass fix, other peaks also match (PR #1440)
        fabs(it->first - 474.248) < 0.001 ||
        fabs(it->first - 545.285) < 0.001)

    {
      TEST_EQUAL(it->second.score > 1, true)
    }
    else
    {
      TEST_EQUAL(it->second.score <= 1, true)
    }
  }


}
END_SECTION

START_SECTION((CompNovoIonScoring& operator=(const CompNovoIonScoring &source)))
{
  CompNovoIonScoring cnis;
  Param p(cnis.getParameters());
  p.setValue("fragment_mass_tolerance", 0.5);
  cnis.setParameters(p);
	CompNovoIonScoring cnis2;
	cnis2 = cnis;
  TEST_EQUAL(cnis2.getParameters() == p, true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



