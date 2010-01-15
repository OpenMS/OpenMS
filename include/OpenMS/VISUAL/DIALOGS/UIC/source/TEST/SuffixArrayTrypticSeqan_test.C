// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <vector>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticSeqan.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SuffixArrayTrypticSeqan, "$Id: SuffixArrayTrypticSeqan_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SuffixArrayTrypticSeqan* ptr = 0;
const String text = "$AAARAA$ARARP$";

SuffixArrayTrypticSeqan* sa = new SuffixArrayTrypticSeqan(text,"");

START_SECTION(SuffixArrayTrypticSeqan(const String &st, const String &filename, const WeightWrapper::WEIGHTMODE weight_mode=WeightWrapper::MONO))
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArrayTrypticSeqan("A",""));
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArrayTrypticSeqan("$A",""));
	ptr = new SuffixArrayTrypticSeqan("$","");
	TEST_NOT_EQUAL(ptr, 0);
	TEST_EXCEPTION (Exception::FileNotFound,new SuffixArrayTrypticSeqan(text,"FileThatNotExists"));
	
END_SECTION

START_SECTION(bool isDigestingEnd(const char aa1, const char aa2) const)
	TEST_EQUAL (sa->isDigestingEnd('R','R'),true);
	TEST_EQUAL (sa->isDigestingEnd('K','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','P'),false);
	TEST_EQUAL (sa->isDigestingEnd('K','P'),false);
	TEST_EQUAL (sa->isDigestingEnd('A','R'),false);
END_SECTION



START_SECTION([EXTRA]SuffixArrayTrypticSeqan::findSpec(const std::vector<DoubleReal> & spec ))
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
	sa = new SuffixArrayTrypticSeqan(text,"");
	vector<DoubleReal> spec;
	//spec.push_back(178.1864 + 18.0);
	//spec.push_back(441.4806 + 18.0);
	//spec.push_back(245.2816);               // AR
	//spec.push_back(387.4392);


	spec.push_back(AASequence("AR").getMonoWeight(Residue::Full));		// AR
	spec.push_back(AASequence("AAAR").getMonoWeight(Residue::Full));		// AAAR

	cerr << 245.2816 << " " << AASequence("AR").getMonoWeight(Residue::Full) << endl;
	cerr << 387.4392 << " " << AASequence("AAAR").getMonoWeight(Residue::Full) << endl;

	vector <vector< pair<pair<SignedSize, SignedSize>,DoubleReal> > > res;
	cerr << "res.size()=" << res.size() << endl;
	sa->findSpec(res, spec);
	TEST_EQUAL(res.size(),spec.size());
	for (Size i = 0; i<res.size();i++)
	{
		TEST_EQUAL(res.at(i).size(),1);
	}
	
	TEST_EQUAL(res.at(0).at(0).first.first,8)
	TEST_EQUAL(res.at(0).at(0).first.second,2)
	TEST_EQUAL(res.at(1).at(0).first.first,1)
	TEST_EQUAL(res.at(1).at(0).first.second,4)
	spec.clear();
	const vector<DoubleReal> specc2 (spec);
	res.clear();
	sa->findSpec(res, specc2);
	TEST_EQUAL(res.size(),0);
	spec.push_back(441.4806);	
	spec.push_back(178.1864);
	const vector<DoubleReal> specc3 (spec);
	res.clear();
	TEST_EXCEPTION(Exception::InvalidValue, sa->findSpec(res, specc3));
	ifstream i_stream;
	i_stream.open(OPENMS_GET_TEST_DATA_PATH("SuffixArrayTrypticSeqan_test.txt"));
	String txt;
	getline(i_stream,txt);
	sa = new SuffixArrayTrypticSeqan(txt,"");
	vector<DoubleReal> spec_new;
	for (int i = 500; i < 5000; i += 197) 
	{
		spec_new.push_back((DoubleReal)i);
	}

	vector<DoubleReal> specc_new(spec_new);
	res.clear();
	sa->findSpec(res, specc_new);
	
	//checking for doubled results;
	for (Size i = 0; i < res.size();++i)
	{
		for (Size j = 0;j<res.at(i).size();++j)
		{
			for (Size k = j+1; k < res.at(i).size();++k)
			{
				TEST_EQUAL(res.at(i).at(j).first.first==res.at(i).at(k).first.first && res.at(i).at(j).first.second==res.at(i).at(k).first.second, false);			
			}
		}
	}

	TOLERANCE_ABSOLUTE(0.55)
	sa->setTolerance(0.5);
	// checking if the mass of the found candidates is correct
	// checking if the next character is not a P
	
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			if (txt[res.at(i).at(j).first.first-1]!='$') TEST_NOT_EQUAL(seq[0],'P');
			if (txt[res.at(i).at(j).first.first+res.at(i).at(j).first.second]!='$') TEST_EQUAL(seq[seq.length()-1]=='R'||seq[seq.length()-1]=='K',true)
			DoubleReal m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length();k++)
			{
				m += masse[(int)seq[k]];
			}
			TEST_REAL_SIMILAR(m,specc_new.at(i));
		}
	}
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
	//cout << "number_of_tags_:" << number_of_tags << endl;
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	vector<String> tags_c(tags);
	sa->setTags(tags_c);
	res.clear();
	sa->findSpec(res, specc_new);
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
			if (!has_tag) cout << seq << endl;
			TEST_EQUAL(has_tag,true);
			TEST_EQUAL(res.at(i).at(j).second, 0);
			
			res_with_tags.push_back(seq);
		}
	}

	for (Size i = 0; i < res_with_tags_exp.size(); ++i)
	{
		bool was_found = false;
		for (Size j = 0; j < res_with_tags.size(); ++j)
		{
			if (res_with_tags_exp.at(i)==res_with_tags.at(j))
			{
				was_found=true;
				break;
			}
		}
		if (!was_found) cout << res_with_tags_exp.at(i) << endl;
	}
	//cout << "mod: 1" << endl;

	sa->setNumberOfModifications(1);
	sa->setUseTags(false);
	res.clear();
	sa->findSpec(res, specc_new);
	
	for (Size i = 0; i < res.size(); i++)
	{
		for (Size j = 0; j < res[i].size(); j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			DoubleReal m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length(); k++)
			{
				m += masse[(int)seq[k]];
			}
			
			TEST_REAL_SIMILAR(m + res[i][j].second, specc_new[i]);
		}
	}
	
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
