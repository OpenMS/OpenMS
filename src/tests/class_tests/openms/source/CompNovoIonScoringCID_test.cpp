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
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoringCID, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIonScoringCID* ptr = nullptr;
CompNovoIonScoringCID* nullPointer = nullptr;
START_SECTION(CompNovoIonScoringCID())
{
 ptr = new CompNovoIonScoringCID();
 TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((CompNovoIonScoringCID(const CompNovoIonScoringCID &source)))
{
  CompNovoIonScoringCID cnis;
  Param p(cnis.getParameters());
  p.setValue("fragment_mass_tolerance", 0.6);
  cnis.setParameters(p);
  TEST_EQUAL(CompNovoIonScoringCID(cnis).getParameters() == p, true)
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoringCID()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void scoreSpectrum(Map<double, IonScore>& CID_ion_scores, PeakSpectrum& CID_spec, double precursor_weight, Size charge)))
{
  Map<double, CompNovoIonScoringBase::IonScore> ion_scores;
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

  CompNovoIonScoringCID cnis;
  cnis.scoreSpectrum(ion_scores, spec, 1018.48, 1);
  
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
    if (fabs(it->first - 903.468292000971) < 0.01 || 
        fabs(it->first - 756.399878084171) < 0.01 ||
        fabs(it->first - 659.347114231171) < 0.01 ||
        fabs(it->first - 546.263050250571) < 0.01 ||
        fabs(it->first - 475.225936461371) < 0.01 ||
        fabs(it->first - 361.183009010571) < 0.01 ||
        fabs(it->first - 304.161545285171) < 0.01 ||
        fabs(it->first - 175.118952187571) < 0.01 ||
        fabs(it->first - 263.102633417371) < 0.01 ||
        fabs(it->first - 360.155397270371) < 0.01 ||
        /*fabs(it->first - 473.239461250971) < 0.01 ||*/
        /*fabs(it->first - 544.276575040171) < 0.01 ||*/
        fabs(it->first - 658.319502490971) < 0.01 ||
        fabs(it->first - 715.340966216371) < 0.01 ||
        fabs(it->first - 844.383559313971) < 0.01)
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

START_SECTION((CompNovoIonScoringCID& operator=(const CompNovoIonScoringCID &source)))
{
  CompNovoIonScoringCID cnis;
  Param p(cnis.getParameters());
  p.setValue("fragment_mass_tolerance", 0.6);
  cnis.setParameters(p);
  CompNovoIonScoringCID cnis2;
  cnis2 = cnis;
  TEST_EQUAL(cnis2.getParameters() == cnis.getParameters(), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

