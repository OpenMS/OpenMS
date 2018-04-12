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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/ExperimentalDesign.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExperimentalDesign, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExperimentalDesign* ptr = 0;
ExperimentalDesign* null_ptr = 0;

ExperimentalDesignFile design = ExperimentalDesignFile::load(OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_1.tsv", false));

START_SECTION(ExperimentalDesign())
{
	ptr = new ExperimentalDesign();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ExperimentalDesign())
{
	delete ptr;
}
END_SECTION

START_SECTION((ExperimentalDesign(MSFileSection msfile_section, SampleSection sample_section)))
{
  // TODO
}
END_SECTION

START_SECTION((const MSFileSection& getMSFileSection() const ))
{
  MSFileSection fs = design.getMSFileSection();
}
END_SECTION

START_SECTION((void setMSFileSection(const MSFileSection &msfile_section)))
{
  ExperimentalDesign design2 = design;
  MSFileSection fs;
  design2.setMSFileSection(fs);
}
END_SECTION

START_SECTION((const SampleSection& getSampleSection() const ))
{
}
END_SECTION

START_SECTION((void setSampleSection(const SampleSection &sample_section)))
{
  ExperimentalDesign design2 = design;
  SampleSection fs;
  design2.setSampleSection(fs);
}
END_SECTION

START_SECTION((std::vector< String > getFileNames(bool basename) const ))
{
  std::vector< String > fns = design.getFileNames();
}
END_SECTION

START_SECTION((std::vector<unsigned> getLabels() const ))
{
  std::vector< unsigned > ls = design.getLabels();
}
END_SECTION

START_SECTION((std::vector<unsigned> getFractions() const ))
{
  std::vector< unsigned > fs = design.getFractions();
}
END_SECTION

START_SECTION((std::map<unsigned int, std::vector<String> > getFractionToMSFilesMapping() const ))
{
  std::map<unsigned int, std::vector<String> > f2ms = design.getFractionToMSFilesMapping();
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToSampleMapping(bool) const ))
{
  std::map< std::pair< String, unsigned > pl2s = getPathLabelToSampleMapping(true);
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionMapping(bool) const ))
{
  std::map< std::pair< String, unsigned > pl2f = getPathLabelToFractionMappingMapping(true);
}
END_SECTION

START_SECTION((std::map< std::pair< String, unsigned >, unsigned> getPathLabelToFractionGroupMapping(bool) const ))
{
  std::map< std::pair< String, unsigned > pl2fg = design.getPathLabelToFractionGroupMapping(true);
}
END_SECTION

START_SECTION((unsigned getNumberOfSamples() const ))
{
  unsigned ns = design.getNumberOfSamples();
}
END_SECTION

START_SECTION((unsigned getNumberOfFractions() const ))
{
  unsigned nf = design.getNumberOfFractions();
}
END_SECTION

START_SECTION((unsigned getNumberOfLabels() const ))
{
  unsigned nl = design.getNumberOfLabels();
}
END_SECTION

START_SECTION((unsigned getNumberOfMSFiles() const ))
{
  unsigned nms = design.getNumberOfMSFiles();
}
END_SECTION

START_SECTION((unsigned getNumberOfFractionGroups() const ))
{
  unsigned nfg = design.getNumberOfFractionGroups(); 
}
END_SECTION

START_SECTION((unsigned getSample(unsigned fraction_group, unsigned label=1)))
{
  unsigned s = design.getSample(1, 1); 
}
END_SECTION

START_SECTION((bool isFractionated() const ))
{
  bool b = design.isFractionated(); 
}
END_SECTION

START_SECTION((bool sameNrOfMSFilesPerFraction() const ))
{
  bool b = design.sameNrOfMSFilesPerFraction(); 
}
END_SECTION

START_SECTION((static ExperimentalDesign fromConsensusMap(const ConsensusMap &c)))
{
  // TODO
}
END_SECTION

START_SECTION((static ExperimentalDesign fromFeatureMap(const FeatureMap &f)))
{
  // TODO
}
END_SECTION

START_SECTION((static ExperimentalDesign fromIdentifications(const std::vector< ProteinIdentification > &proteins)))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] SampleSection(std::vector< std::vector< String > > _content, std::map< unsigned, Size > _sample_to_rowindex, std::map< String, Size > _columnname_to_columnindex)))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] std::set< unsigned > getSamples() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] std::set< String > getFactors() const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] bool hasSample(unsigned sample) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] bool hasFactor(const String &factor) const ))
{
  // TODO
}
END_SECTION

START_SECTION(([ExperimentalDesign::SampleSection] String getFactorValue(unsigned sample, const String &factor)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



