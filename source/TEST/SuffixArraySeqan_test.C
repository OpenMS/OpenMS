// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
#include <vector>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArraySeqan.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SuffixArraySeqan, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SuffixArraySeqan* ptr = 0;
const String text = "$AAARAA$ARARP$";

SuffixArraySeqan* sa = new SuffixArraySeqan(text,"");

START_SECTION((SuffixArraySeqan(const String &st, const String &filename, const WeightWrapper::WEIGHTMODE weight_mode=WeightWrapper::MONO)))
{
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArraySeqan("A",""));
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArraySeqan("$A",""));
	ptr = new SuffixArraySeqan("$","");
	TEST_NOT_EQUAL(ptr, 0);
	TEST_EXCEPTION (Exception::FileNotFound,new SuffixArraySeqan(text,"FileThatNotExists"));
}
END_SECTION

START_SECTION((SuffixArraySeqan(const SuffixArraySeqan &source)))
{
  SuffixArraySeqan* sa_new = new SuffixArraySeqan(text,"");
	sa_new->setTolerance(0.1);
	sa_new->setNumberOfModifications(1);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	sa_new->setTags(tags);
	SuffixArraySeqan sa2 (*sa_new);
	TEST_EQUAL (sa_new->getTolerance(),sa2.getTolerance());
	TEST_EQUAL (sa_new->getNumberOfModifications(),sa2.getNumberOfModifications());
	TEST_EQUAL (sa_new->getUseTags(),sa2.getUseTags());
	TEST_EQUAL (sa_new->getTags().size(),sa2.getTags().size());
	for (Size i = 0; i < sa2.getTags().size();++i)
	{
		TEST_EQUAL (sa_new->getTags().at(i),sa2.getTags().at(i));
	}
}
END_SECTION

START_SECTION((~SuffixArraySeqan()))
{
  delete ptr;
}
END_SECTION

START_SECTION(void printStatistic())
{
	NOT_TESTABLE;
}
END_SECTION

START_SECTION(bool isDigestingEnd(const char aa1, const char aa2) const )
{
	TEST_EQUAL (sa->isDigestingEnd('R','R'),true);
	TEST_EQUAL (sa->isDigestingEnd('K','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','P'),true);
	TEST_EQUAL (sa->isDigestingEnd('K','P'),true);
	TEST_EQUAL (sa->isDigestingEnd('A','R'),true);
}
END_SECTION

START_SECTION(DoubleReal getTolerance() const )
{
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
}
END_SECTION

START_SECTION(void setTolerance(DoubleReal t))
{
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
}
END_SECTION

START_SECTION(Size getNumberOfModifications ())
{
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
}
END_SECTION

START_SECTION(String toString())
{
	SuffixArraySeqan new_sa(text, "");
	String sa_string = new_sa.toString();
	// not implemented in this SA, hence string is empty
	TEST_STRING_EQUAL(sa_string, "");
}
END_SECTION

START_SECTION(void setNumberOfModifications(Size number_of_mods))
{
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
	TEST_EXCEPTION(Exception::InvalidValue,sa->setTolerance(-0.5));
}
END_SECTION

START_SECTION(void setTags(const std::vector< OpenMS::String > &tags))
{
	SuffixArraySeqan * satc = new SuffixArraySeqan(text,"");
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	satc->setTags(tags);
	vector<String> res = satc->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
}
END_SECTION

START_SECTION(const std::vector<OpenMS::String>& getTags())
{
	SuffixArraySeqan * satc = new SuffixArraySeqan(text,"");
	TEST_EQUAL(satc->getTags().size(),0);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	vector<String> res = satc->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
}
END_SECTION

START_SECTION((void setUseTags(bool use_tags)))
{
	SuffixArraySeqan * satc = new SuffixArraySeqan(text,"");
	TEST_EQUAL(satc->getUseTags(),false);
	satc->setUseTags(1);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	satc->setUseTags(0);
	TEST_EQUAL(satc->getUseTags(),false);
}
END_SECTION

START_SECTION(bool getUseTags())
{
	SuffixArraySeqan * satc = new SuffixArraySeqan(text,"");
	TEST_EQUAL(satc->getUseTags(),false);
	satc->setUseTags(1);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	satc->setUseTags(0);
	TEST_EQUAL(satc->getUseTags(),false);
}
END_SECTION

START_SECTION((bool open(const String &filename)))
{
	TEST_EXCEPTION (Exception::FileNotFound,sa->open("FileThatNotExists"));
	NOT_TESTABLE // will be tested in next test
}
END_SECTION


START_SECTION((bool save(const String &filename)))
{
	String filename;
	NEW_TMP_FILE(filename)
	sa->save(filename);
	SuffixArraySeqan sa2("$", filename);
	sa2.open(filename);
	NOT_TESTABLE
}
END_SECTION

#if 1

START_SECTION((void findSpec(std::vector< std::vector< std::pair< std::pair< SignedSize, SignedSize >, DoubleReal > > > &candidates, const std::vector< DoubleReal > &spec)))
{
	DoubleReal masse[255];
	ResidueDB* rdb = ResidueDB::getInstance();
		
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
		
	for (Size i = 0; i<255;++i)
	{
		masse[i]=0;
	}
	for (Size i = 0; i<strlen(aa);++i)
	{
		const Residue * r = rdb->getResidue(aa[i]);
		masse[(int)aa[i]]=r->getMonoWeight(Residue::Internal);
	}
	sa = new SuffixArraySeqan(text,"");
	vector<DoubleReal> spec;
	//spec.push_back(178.1864);
	//spec.push_back(441.4806);
	spec.push_back(AASequence("AR").getMonoWeight(Residue::Full));
  spec.push_back(AASequence("AAAR").getMonoWeight(Residue::Full));
	const vector<DoubleReal> specc (spec);
	vector <vector< pair<pair<SignedSize, SignedSize>,DoubleReal> > > res;
	sa->findSpec(res, specc);
	TEST_EQUAL(res.size(),specc.size());
	
	TEST_EQUAL(res[0].size(), 5);
	TEST_EQUAL(res[1].size(), 3);

	TEST_EQUAL(res.at(0).at(0).first.first,8)
	TEST_EQUAL(res.at(0).at(0).first.second,2)
	TEST_EQUAL(res.at(1).at(0).first.first,1)
	TEST_EQUAL(res.at(1).at(0).first.second,4)

	spec.clear();
	vector<DoubleReal> specc2 (spec);
	res.clear();
	sa->findSpec(res, specc2);
	TEST_EQUAL(res.size(),0);
	spec.push_back(441.4806);	
	spec.push_back(178.1864);
	const vector<DoubleReal> specc3 (spec);
	res.clear();
	TEST_EXCEPTION(Exception::InvalidValue, sa->findSpec(res, specc3));
	ifstream i_stream;
	i_stream.open(OPENMS_GET_TEST_DATA_PATH("SuffixArraySeqan_test.txt"));
	String txt;
	getline(i_stream,txt);
	sa = new SuffixArraySeqan(txt,"");
	STATUS("Okay!");
	vector<DoubleReal> spec_new;
	for (int i = 500; i < 5000; i += 197)
	{
		spec_new.push_back((DoubleReal)i);
	}
	STATUS("Okay!");
	vector<DoubleReal> specc_new(spec_new);
	STATUS("Okay!");
	res.clear();
	STATUS("Okay!");
	sa->findSpec(res, specc_new);
	//checking for doubled results;
	STATUS("Okay!");
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			for (Size k = j+1; k < res.at(i).size();k++)
			{
				TEST_EQUAL(res.at(i).at(j).first.first==res.at(i).at(k).first.first && res.at(i).at(j).first.second==res.at(i).at(k).first.second, false);
				
			}
		}
	}
	STATUS("Okay!");
	TOLERANCE_ABSOLUTE(0.55)
	sa->setTolerance(0.5);
	// checking if the mass of the found candidates is correct
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			
			DoubleReal m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length();k++)
			{
				m+=masse[(int)seq[k]];
			}
			TEST_REAL_SIMILAR(m,specc_new.at(i));
		}
	}
	STATUS("Okay!");

	// getting all candidates with tags 
	Size number_of_tags=0;
	vector<String> res_with_tags_exp;
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			bool has_tag = false;
			for (Size k = 2; k < seq.length();k++)
			{
				if (seq.substr(k-2,3)=="AAA"||seq.substr(k-2,3)=="ARA")
				{
					has_tag=true;
					break;
				}
			}
			if (has_tag) {
				++number_of_tags;
				res_with_tags_exp.push_back(seq);
			}
			
		}
	}

	STATUS("Okay!");

	//cout << "number_of_tags_:" << number_of_tags << endl;
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c (tags);
	sa->setTags(tags_c);
	res.clear();
	sa->findSpec(res, specc_new);
	STATUS("Okay!");
	vector<String> res_with_tags;
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			bool has_tag = false;
			for (Size k = 2; k < seq.length();k++)
			{
				if (seq.substr(k-2,3)=="AAA"||seq.substr(k-2,3)=="ARA")
				{
					has_tag=true;
					break;
				}
			}
			if (!has_tag) cout <<seq << endl;
			TEST_EQUAL(has_tag, true);
			TEST_EQUAL(res.at(i).at(j).second,0);
			
			res_with_tags.push_back(seq);
		}
	}
	STATUS("Okay!");
	for (Size i = 0; i < res_with_tags_exp.size();++i){
		bool was_found = false;
		for (Size j = 0; j < res_with_tags.size();++j){
			if (res_with_tags_exp.at(i)==res_with_tags.at(j)){
				was_found=true;
				break;
			}
		}
		if (!was_found) cout << res_with_tags_exp.at(i) << endl;
	}
	// cout << "mod: 1" << endl;
	sa->setNumberOfModifications(1);
	sa->setUseTags(false);
	res.clear();
	sa->findSpec(res, specc_new);
	
	STATUS("Okay!");

	//Checking if mass is correct
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			DoubleReal m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length();k++)
			{
				m+=masse[(int)seq[k]];
			}
			TEST_REAL_SIMILAR(m+res.at(i).at(j).second,specc_new.at(i));
			
		}
	}
	spec.clear();
	STATUS("Okay!");

	//testing if a candidate can belong to serveal input masses
	spec.push_back(441.4806);
	spec.push_back(441.4806);
	const vector<DoubleReal> specc4 (spec);
	sa->setNumberOfModifications(0);
	sa->setUseTags(false);
	res.clear();
	sa->findSpec(res, specc4);
	TEST_EQUAL(res.at(0).size(),res.at(1).size());
	for (Size j = 0; j < res.at(0).size();++j)
	{
		TEST_EQUAL(res.at(0).at(j).first.first,res.at(1).at(j).first.first);
		TEST_EQUAL(res.at(0).at(j).first.second,res.at(1).at(j).first.second);
		TEST_EQUAL(res.at(0).at(j).second,res.at(1).at(j).second);
	}
	STATUS("Okay!");

}
END_SECTION

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
