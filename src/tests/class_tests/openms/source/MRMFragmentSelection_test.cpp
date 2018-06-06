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
#include <OpenMS/ANALYSIS/MRM/MRMFragmentSelection.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMFragmentSelection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMFragmentSelection* ptr = nullptr;
MRMFragmentSelection* nullPointer = nullptr;
START_SECTION(MRMFragmentSelection())
{
  ptr = new MRMFragmentSelection();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MRMFragmentSelection())
{
  delete ptr;
}
END_SECTION

START_SECTION((MRMFragmentSelection(const MRMFragmentSelection &rhs)))
{
  MRMFragmentSelection mrmfs;
  Param p = mrmfs.getParameters();
  p.setValue("num_top_peaks", 18);
  mrmfs.setParameters(p);
  TEST_EQUAL(MRMFragmentSelection(mrmfs).getParameters() == p, true)
}
END_SECTION

START_SECTION((MRMFragmentSelection& operator=(const MRMFragmentSelection &rhs)))
{
  MRMFragmentSelection mrmfs;
  Param p = mrmfs.getParameters();
  p.setValue("num_top_peaks", 18);
  mrmfs.setParameters(p);
  MRMFragmentSelection mrmfs2;
  mrmfs2 = mrmfs;
  TEST_EQUAL(mrmfs2.getParameters() == p, true)
}
END_SECTION

START_SECTION((void selectFragments(std::vector< Peak1D > &selected_peaks, const PeakSpectrum &spec)))
{
  PeakSpectrum spec;
  TheoreticalSpectrumGenerator tsg;
  Param tsg_param(tsg.getParameters());
  tsg_param.setValue("add_metainfo", "true");
  tsg.setParameters(tsg_param);
  tsg.getSpectrum(spec, AASequence::fromString("DFPIANGER"), 1, 1);

  spec.sortByPosition();
  Precursor prec;
  prec.setMZ(1019.1);
  vector<Precursor> precursors;
  precursors.push_back(prec);
  spec.setPrecursors(precursors);

  PeptideHit hit;
  hit.setCharge(1);
  hit.setSequence(AASequence::fromString("DFPIANGER"));
  vector<PeptideHit> hits;
  hits.push_back(hit);
  PeptideIdentification id;
  id.setHits(hits);
  vector<PeptideIdentification> ids;
  ids.push_back(id);
  spec.setPeptideIdentifications(ids);

  MRMFragmentSelection mrmfs;
  Param p(mrmfs.getParameters());
  p.setValue("num_top_peaks", 1);
  p.setValue("allowed_ion_types", ListUtils::create<String>("y"));
  mrmfs.setParameters(p);

  vector<Peak1D> selected_peaks;
  mrmfs.selectFragments(selected_peaks, spec);
  TEST_EQUAL(selected_peaks.size(), 1)

  p.setValue("num_top_peaks", 3);
  p.setValue("min_pos_precursor_percentage", 10.0);
  mrmfs.setParameters(p);
  selected_peaks.clear();
  mrmfs.selectFragments(selected_peaks, spec);
  TEST_EQUAL(selected_peaks.size(), 3)

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



