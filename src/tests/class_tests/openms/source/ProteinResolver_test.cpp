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
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProteinResolver, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ProteinResolver* ptr = nullptr;
ProteinResolver* null_ptr = nullptr;
START_SECTION(ProteinResolver())
{
  ptr = new ProteinResolver();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ProteinResolver())
{
  delete ptr;
}
END_SECTION

START_SECTION((ProteinResolver(const ProteinResolver &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~ProteinResolver()))
{
  // TODO
}
END_SECTION

START_SECTION((ProteinResolver& operator=(const ProteinResolver &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void resolve(std::vector< ISDGroup > &isd_groups, std::vector< MSDGroup > &msd_groups, std::vector< FASTAFile::FASTAEntry > &protein_data, std::vector< Size > &reindexed_proteins, std::vector< Size > &reindexed_peptides, std::vector< ProteinEntry > &protein_nodes, std::vector< PeptideEntry > &peptide_nodes, std::vector< PeptideIdentification > peptide_identifications, ConsensusMap consensus, EnzymaticDigestion digestor, bool id, UInt min_size)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinsAndPeptidesmzTab(std::vector< ProteinEntry > &protein_nodes, std::vector< PeptideEntry > &peptide_nodes, std::vector< Size > &reindexed_proteins, std::vector< Size > &reindexed_peptides, std::vector< PeptideIdentification > &peptide_identifications)))
{
  // TODO
}
END_SECTION

START_SECTION((void writePeptideTable(std::vector< PeptideEntry > &peptides, std::vector< Size > &reindexed_peptides, std::vector< PeptideIdentification > &identifications, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writePeptideTable(std::vector< PeptideEntry > &peptides, std::vector< Size > &reindexed_peptides, ConsensusMap &consensus, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinTable(std::vector< ProteinEntry > &proteins, std::vector< Size > &reindexed_proteins, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void writeProteinGroups(std::vector< ISDGroup > &isd_groups, std::vector< MSDGroup > &msd_groups, String &output_file)))
{
  // TODO
}
END_SECTION

START_SECTION((void countTargetDecoy(std::vector< MSDGroup > &msd_groups, ConsensusMap &consensus)))
{
  // TODO
}
END_SECTION

START_SECTION((void countTargetDecoy(std::vector< MSDGroup > &msd_groups, std::vector< PeptideIdentification > &peptide_nodes)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
