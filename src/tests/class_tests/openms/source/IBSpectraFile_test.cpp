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

#include <OpenMS/FORMAT/IBSpectraFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

START_TEST(IBSpectraFile, "$Id$")

IBSpectraFile* ptr = nullptr;
IBSpectraFile* nullPointer = nullptr;

START_SECTION((IBSpectraFile()))
{
  ptr = new IBSpectraFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((IBSpectraFile(const IBSpectraFile& other)))
{
  IBSpectraFile ibfile(*ptr);
  TEST_NOT_EQUAL(&ibfile, nullPointer)
}
END_SECTION

START_SECTION((IBSpectraFile& operator=(const IBSpectraFile& rhs)))
{
  IBSpectraFile ibfile;
  ibfile = *ptr;
  TEST_NOT_EQUAL(&ibfile, nullPointer)
}
END_SECTION

START_SECTION((void store(const String& filename, const ConsensusMap& cm)))
{
  // test invalid ConsensusMap
  ConsensusMap cm_no_ms2quant;
  cm_no_ms2quant.setExperimentType("not-isobaric");

  IBSpectraFile ibfile_no_ms2quant;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ibfile_no_ms2quant.store("not-a-file-name", cm_no_ms2quant), "Given ConsensusMap does not hold any isobaric quantification data.")

  // test wrong channel count
  ConsensusMap cm_wrong_channel_count;
  cm_wrong_channel_count.setExperimentType("labeled_MS2");
  ConsensusMap::FileDescription channel1;
  ConsensusMap::FileDescription channel2;
  ConsensusMap::FileDescription channel3;
  cm_wrong_channel_count.getFileDescriptions()[0] = channel1;
  cm_wrong_channel_count.getFileDescriptions()[1] = channel2;
  cm_wrong_channel_count.getFileDescriptions()[2] = channel3;

  IBSpectraFile ibfile_wrong_channel_count;
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, ibfile_wrong_channel_count.store("not-a-file-name", cm_wrong_channel_count), "Could not guess isobaric quantification data from ConsensusMap due to non-matching number of input maps.")

  // test a real example
  ConsensusMap cm;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("IBSpectraFile.consensusXML"),cm);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);

  IBSpectraFile ibfile;
  ibfile.store(tmp_filename, cm);

  TEST_FILE_SIMILAR(tmp_filename.c_str(), OPENMS_GET_TEST_DATA_PATH("IBSpectraFile.ibspectra.csv"))
}
END_SECTION

END_TEST
