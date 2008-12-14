// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl,Andreas Bertsch$
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


START_SECTION(SuffixArrayPeptideFinder(const String& filename, const String& method))
	ptr = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","seqan");
	ptr = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	ptr = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticCompressed");
	TEST_EXCEPTION(Exception::InvalidValue,new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","bla"));
	TEST_EXCEPTION(Exception::FileNotFound,new SuffixArrayPeptideFinder("FileThatNotExists","seqan"));
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION


START_SECTION(~SuffixArrayPeptideFinder())
	delete ptr;
END_SECTION

START_SECTION(SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder &source))
	ptr = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","seqan");
	ptr->setTolerance(0.1);
	ptr->setNumberOfModifications(2);
	SuffixArrayPeptideFinder * new_ptr = new SuffixArrayPeptideFinder (*ptr);
	TEST_EQUAL (ptr->getNumberOfModifications(),new_ptr->getNumberOfModifications());
	TEST_EQUAL (ptr->getTolerance(),new_ptr->getTolerance());
END_SECTION

START_SECTION(double getTolerance() const)
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
END_SECTION

START_SECTION(void setTolerance(const double t))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
	TEST_EXCEPTION(Exception::InvalidValue,sa->setTolerance(-0.5));
END_SECTION

START_SECTION(void setNumberOfModifications(UInt number_of_mods) const)
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
END_SECTION

START_SECTION(UInt getNumberOfModifications() const)
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
END_SECTION

START_SECTION(void setTags(const std::vector< OpenMS::String > &tags))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa->setTags(tags);
	std::vector<String> res = sa->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
END_SECTION

START_SECTION(const std::vector<OpenMS::String>& getTags())
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_EQUAL(sa->getTags().size(),0);
	TEST_EQUAL(sa->getUseTags(),0);
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa->setTags(tags);
	TEST_EQUAL(sa->getUseTags(),1);
	std::vector<String> res = sa->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
END_SECTION

START_SECTION(void setUseTags(bool use_tags))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_EQUAL(sa->getUseTags(),0);
	sa->setUseTags(1);
	TEST_EQUAL(sa->getUseTags(),0);
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa->setTags(tags);
	TEST_EQUAL(sa->getUseTags(),1);
	sa->setUseTags(0);
	TEST_EQUAL(sa->getUseTags(),0);
END_SECTION

START_SECTION(bool getUseTags())
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	TEST_EQUAL(sa->getUseTags(),0);
	sa->setUseTags(1);
	TEST_EQUAL(sa->getUseTags(),0);
	std::vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const std::vector<String> tags_c (tags);
	sa->setTags(tags);
	TEST_EQUAL(sa->getUseTags(),1);
	sa->setUseTags(0);
	TEST_EQUAL(sa->getUseTags(),0);
END_SECTION

START_SECTION(void setModificationOutputMethod(const String &s))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","seqan");
	sa->setModificationOutputMethod ("stringChecked");
	TEST_EQUAL("stringChecked",sa->getModificationOutputMethod());
	sa->setModificationOutputMethod ("stringUnchecked");
	TEST_EQUAL("stringUnchecked",sa->getModificationOutputMethod());
	sa->setModificationOutputMethod ("mass");
	TEST_EQUAL("mass",sa->getModificationOutputMethod());
	TEST_EXCEPTION(Exception::InvalidValue,sa->setModificationOutputMethod ("bla"));
END_SECTION

START_SECTION(String getModificationOutputMethod())
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","seqan");
	TEST_EQUAL("mass",sa->getModificationOutputMethod());
	sa->setModificationOutputMethod ("stringChecked");
	TEST_EQUAL("stringChecked",sa->getModificationOutputMethod());
	sa->setModificationOutputMethod ("stringUnchecked");
	TEST_EQUAL("stringUnchecked",sa->getModificationOutputMethod());
	sa->setModificationOutputMethod ("mass");
	TEST_EQUAL("mass",sa->getModificationOutputMethod());
END_SECTION

START_SECTION((void getCandidates(std::vector< std::vector< std::pair< FASTAEntry, String > > > &candidates, const std::vector< double > &spec)))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	vector<double> spec;
	spec.push_back(178.1864);
	spec.push_back(441.4806);
	const vector<double> specc (spec);
	sa->setTolerance(0.5);
	sa->setNumberOfModifications(0);
	vector<vector<pair<FASTAEntry,String> > > res2;
	sa->getCandidates(res2, specc);
	/*
	for (vector<vector<pair<FASTAEntry,String> > >::const_iterator it1 = res2.begin(); it1 != res2.end(); ++it1)
	{
		for (vector<pair<FASTAEntry, String> >::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2)
		{
			cerr << it2->first.first << " ##### " << it2->first.second << " " << AASequence(it2->first.second).getAverageWeight() << endl;
		}
	}
	*/
	FastaIterator * fit = new FastaIterator ();
	fit->setFastaFile("data/SuffixArrayPeptideFinder_test.fasta");
	fit->begin();
	std::map<String,String> fasta_map;
	while (!fit->isAtEnd())
	{
		fasta_map[(**fit).first]=(**fit).second;
		++*fit;
	}
	for (unsigned int i = 0; i < res2.size(); ++i)
	{
		for (unsigned int j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			unsigned int l = pep_seq.length();
			bool found = false;
			for (unsigned int k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,1);
			//if (!found ) std::cout<<pep_seq <<":"<<complete_seq<<std::endl;
			TEST_EQUAL (res2.at(i).at(j).second,"");
		}
	}
	sa->setNumberOfModifications(1);
	res2.clear();
	sa->getCandidates(res2, specc);
	for (unsigned int i = 0; i < res2.size(); ++i)
	{
		for (unsigned int j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			unsigned int l = pep_seq.length();
			bool found = false;
			for (unsigned int k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,1);
			String mod_mass = res2.at(i).at(j).second;
			TEST_EQUAL(mod_mass==""||mod_mass=="-1.00794"||mod_mass=="59.044"||mod_mass=="80.9878"||mod_mass=="15.9994",1)
		}
	}
	sa->setModificationOutputMethod ("stringChecked");
	res2.clear();
	sa->getCandidates(res2, specc);
	for (unsigned int i = 0; i < res2.size(); ++i)
	{
		for (unsigned int j = 0; j < res2.at(i).size();++j)
		{
			String pep_seq = res2.at(i).at(j).first.second;
			String complete_seq = fasta_map[res2.at(i).at(j).first.first];
			unsigned int l = pep_seq.length();
			bool found = false;
			for (unsigned int k = l;k<=complete_seq.length();++k)
			{
				found |= complete_seq.substr(k-l,l)==pep_seq;
			}
			TEST_EQUAL (found,1);
			String mod_mass = res2.at(i).at(j).second;
			TEST_EQUAL(mod_mass==""||mod_mass=="[C]"||mod_mass=="[S]"||mod_mass=="[Y]",1)
		}
	}
END_SECTION

START_SECTION(void getCandidates(std::vector< std::vector< std::pair< FASTAEntry, String > > > &candidates, const String &DTA_file))
	SuffixArrayPeptideFinder* sa = new SuffixArrayPeptideFinder("data/SuffixArrayPeptideFinder_test.fasta","trypticSeqan");
	vector<vector<pair<FASTAEntry, String> > > candidates;
	sa->getCandidates(candidates, "data/DTAFile_test.dta");
	TEST_EQUAL(candidates.size(), 25)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



