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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FastaIterator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArrayPeptideFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef std::pair <String, String> FASTAEntry;

START_TEST(SuffixArrayPeptideFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SuffixArrayPeptideFinder* ptr = 0;
SuffixArrayPeptideFinder* nullPointer = 0;


START_SECTION(SuffixArrayPeptideFinder(const String& filename, const String& method, const WeightWrapper::WEIGHTMODE weight_mode=WeightWrapper::MONO))
	ptr = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"seqan");
	ptr = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticSeqan");
	ptr = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticCompressed");
	TEST_EXCEPTION(Exception::InvalidValue,new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"bla"));
	TEST_EXCEPTION(Exception::FileNotFound,new SuffixArrayPeptideFinder("FileThatNotExists","seqan"));
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


SuffixArrayPeptideFinder* sa_tryptic_seqan = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticSeqan");
SuffixArrayPeptideFinder* sa_seqan = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"seqan");


START_SECTION(~SuffixArrayPeptideFinder())
	delete ptr;
END_SECTION

START_SECTION(SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder &source))
	TEST_REAL_SIMILAR (sa_seqan->getTolerance(),0.5);
	TEST_EQUAL (sa_seqan->getNumberOfModifications(),0);
	sa_seqan->setTolerance(0.1);
	sa_seqan->setNumberOfModifications(2);
	SuffixArrayPeptideFinder * new_ptr = new SuffixArrayPeptideFinder (*sa_seqan);
	TEST_EQUAL (sa_seqan->getNumberOfModifications(),new_ptr->getNumberOfModifications());
	TEST_EQUAL (sa_seqan->getTolerance(),new_ptr->getTolerance());
END_SECTION

START_SECTION(DoubleReal getTolerance() const)
	sa_seqan->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa_seqan->getTolerance(),0.1);
	sa_seqan->setTolerance(0.5);
	TEST_REAL_SIMILAR (sa_seqan->getTolerance(),0.5);
END_SECTION

START_SECTION(void setTolerance(const DoubleReal t))
	TEST_EXCEPTION(Exception::InvalidValue,sa_seqan->setTolerance(-0.5));
END_SECTION

START_SECTION(void setNumberOfModifications(Size number_of_mods) const)
	sa_tryptic_seqan->setNumberOfModifications(1);
	TEST_EQUAL (sa_tryptic_seqan->getNumberOfModifications(),1);
	sa_tryptic_seqan->setNumberOfModifications(0);
	TEST_EQUAL (sa_tryptic_seqan->getNumberOfModifications(),0);
END_SECTION

START_SECTION(Size getNumberOfModifications() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setTags(const std::vector< OpenMS::String > &tags))
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(const std::vector<OpenMS::String>& getTags())
	TEST_EQUAL(sa_tryptic_seqan->getTags().size(), 0);
	TEST_EQUAL(sa_tryptic_seqan->getUseTags(), false);
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa_tryptic_seqan->setTags(tags);
	TEST_EQUAL(sa_tryptic_seqan->getUseTags(), true);
	std::vector<String> res = sa_tryptic_seqan->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
END_SECTION

START_SECTION(void setUseTags(bool use_tags))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticSeqan");
	TEST_EQUAL(sa->getUseTags(),false);
	sa->setUseTags(1);
	TEST_EQUAL(sa->getUseTags(),false);
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa->setTags(tags);
	TEST_EQUAL(sa->getUseTags(),true);
	sa->setUseTags(0);
	TEST_EQUAL(sa->getUseTags(),false);
END_SECTION

START_SECTION(bool getUseTags())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setModificationOutputMethod(const String &s))
	TEST_EQUAL(sa_seqan->getModificationOutputMethod(),"mass");
	sa_seqan->setModificationOutputMethod ("stringChecked");
	TEST_EQUAL(sa_seqan->getModificationOutputMethod(),"stringChecked");
	sa_seqan->setModificationOutputMethod ("stringUnchecked");
	TEST_EQUAL(sa_seqan->getModificationOutputMethod(),"stringUnchecked");
	sa_seqan->setModificationOutputMethod ("mass");
	TEST_EQUAL(sa_seqan->getModificationOutputMethod(),"mass");
	TEST_EXCEPTION(Exception::InvalidValue,sa_seqan->setModificationOutputMethod ("bla"));
END_SECTION

START_SECTION(String getModificationOutputMethod())
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void getCandidates(std::vector< std::vector< std::pair< FASTAEntry, String > > > &candidates, const std::vector< DoubleReal > &spec)))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticSeqan");
	vector<DoubleReal> spec;
	spec.push_back(178.1864);
	spec.push_back(441.4806);
	const vector<DoubleReal> specc (spec);
	sa->setTolerance(0.5);
	sa->setNumberOfModifications(0);
	vector<vector<pair<FASTAEntry,String> > > res2;
	sa->getCandidates(res2, specc);
	/*
	for (vector<vector<pair<FASTAEntry,String> > >::const_iterator it1 = res2.begin(); it1 != res2.end(); ++it1)
	{
		for (vector<pair<FASTAEntry, String> >::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
		{
			cerr << it2->first.first << " ##### " << it2->first.second << " " << AASequence(it2->first.second).getMonoWeight() << endl;
		}
	}
	*/
	FastaIterator * fit = new FastaIterator ();
	fit->setFastaFile(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"));
	fit->begin();
	std::map<String,String> fasta_map;
	while (!fit->isAtEnd())
	{
		fasta_map[(**fit).first]=(**fit).second;
		++*fit;
	}
	for (Size i = 0; i < res2.size(); ++i)
	{
		for (Size j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			Size l = pep_seq.length();
			bool found = false;
			for (Size k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,true);
			//if (!found ) std::cout<<pep_seq <<":"<<complete_seq<<std::endl;
			TEST_EQUAL (res2.at(i).at(j).second,"");
		}
	}
	sa->setNumberOfModifications(1);
	res2.clear();
	sa->getCandidates(res2, specc);
	for (Size i = 0; i < res2.size(); ++i)
	{
		for (Size j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			Size l = pep_seq.length();
			bool found = false;
			for (Size k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,true);
			String mod_mass = res2.at(i).at(j).second;
			TEST_EQUAL(mod_mass==""||mod_mass=="-1.00794"||mod_mass=="59.044"||mod_mass=="80.9878"||mod_mass=="15.9994",true)
		}
	}
	sa->setModificationOutputMethod ("stringChecked");
	res2.clear();
	sa->getCandidates(res2, specc);
	for (Size i = 0; i < res2.size(); ++i)
	{
		for (Size j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			Size l = pep_seq.length();
			bool found = false;
			for (Size k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,true);
			String mod_mass = res2.at(i).at(j).second;
			TEST_EQUAL(mod_mass==""||mod_mass=="[C]"||mod_mass=="[S]"||mod_mass=="[Y]",true)
		}
	}
END_SECTION

START_SECTION(void getCandidates(std::vector< std::vector< std::pair< FASTAEntry, String > > > &candidates, const String &DTA_file))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder(OPENMS_GET_TEST_DATA_PATH("SuffixArrayPeptideFinder_test.fasta"),"trypticSeqan");
	vector<vector<pair<FASTAEntry, String> > > candidates;
	sa->getCandidates(candidates, OPENMS_GET_TEST_DATA_PATH("DTAFile_test.dta"));
	TEST_EQUAL(candidates.size(), 25)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



