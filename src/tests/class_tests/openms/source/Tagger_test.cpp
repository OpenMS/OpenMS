// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace OpenMS;

START_TEST(Tagger, "$Id$")

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags))

  TheoreticalSpectrumGenerator tsg;
  Param param = tsg.getParameters();
  param.setValue("add_metainfo", "false");
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "true");
  param.setValue("add_losses", "true");
  param.setValue("add_precursor_peaks", "true");
  tsg.setParameters(param);

  AASequence test_sequence = AASequence::fromString("PEPTIDETESTTHISTAGGER");
  PeakSpectrum spec;
  tsg.getSpectrum(spec, test_sequence, 1, 5);
  Tagger tagger = Tagger(2, 10, 5, 1, 5);
  std::set<std::string> tags;

  TEST_EQUAL(spec.size(), 888);
  tagger.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 11002);

  // runtime benchmark, research tags many times in the same spectrum
  // for (int i = 0; i < 30; i++)
  // {
  //   tags.clear();
  //   tagger.getTag(spec, tags);
  // }

  // write out found tags if necessary
  // for (const std::string& tag : tags)
  // {
  //   std::cout << "TEST TAG: " << tag << std::endl;
  // }

END_SECTION

END_TEST
