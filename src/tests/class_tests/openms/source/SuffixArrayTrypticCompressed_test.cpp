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
#include <OpenMS/test_config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
#include <vector>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArrayTrypticCompressed.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SuffixArrayTrypticCompressed, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SuffixArrayTrypticCompressed* ptr = 0;
SuffixArrayTrypticCompressed* nullPointer = 0;
const String text = "$AAARAA$ARARP$";

SuffixArrayTrypticCompressed* sa = new SuffixArrayTrypticCompressed(text, "");

START_SECTION(SuffixArrayTrypticCompressed(const String &st, const String &filename, const WeightWrapper::WEIGHTMODE weight_mode=WeightWrapper::MONO))
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArrayTrypticCompressed("A",""));
	TEST_EXCEPTION (Exception::InvalidValue,new SuffixArrayTrypticCompressed("$A",""));
	ptr = new SuffixArrayTrypticCompressed("$","");
	TEST_EQUAL(ptr->toString(),"lcp: 0\nskip: 0");
	String s = sa->toString();
	String sTree = s.substr(0,19);
	String lcp = s.substr(25, 4);
	String skip = s.substr(36,4);
	TEST_EQUAL(sTree,"AA\nAAARAA\nARARP\nARP");
	TEST_EQUAL(lcp,"2120");
	TEST_EQUAL(skip,"1210");
  TEST_NOT_EQUAL(ptr, nullPointer);
	TEST_EXCEPTION (Exception::FileNotFound,ptr = new SuffixArrayTrypticCompressed(text,"FileThatNotExists"));
END_SECTION

START_SECTION(SuffixArrayTrypticCompressed(const SuffixArrayTrypticCompressed & sa))
        SuffixArrayTrypticCompressed sa2 (*sa);
	TEST_EQUAL (sa->toString(),sa2.toString());
END_SECTION

START_SECTION(~SuffixArrayTrypticCompressed())
        delete ptr;
END_SECTION

START_SECTION(bool isDigestingEnd(const char aa1, const char aa2) const )
	TEST_EQUAL (sa->isDigestingEnd('R','R'),true);
	TEST_EQUAL (sa->isDigestingEnd('K','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','K'),true);
	TEST_EQUAL (sa->isDigestingEnd('R','P'),false);
	TEST_EQUAL (sa->isDigestingEnd('K','P'),false);
	TEST_EQUAL (sa->isDigestingEnd('A','R'),false);
END_SECTION

START_SECTION(double getTolerance () const)
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
END_SECTION

START_SECTION(void setTolerance(double t))
	TEST_REAL_SIMILAR (sa->getTolerance(),0.5);
	sa->setTolerance(0.1);
	TEST_REAL_SIMILAR (sa->getTolerance(),0.1);
	sa->setTolerance(0.5);
	TEST_EXCEPTION(Exception::InvalidValue,sa->setTolerance(-0.5));
END_SECTION

START_SECTION(Size getNumberOfModifications())
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
END_SECTION

START_SECTION(void setNumberOfModifications(Size number_of_mods))
	TEST_EQUAL (sa->getNumberOfModifications(),0);
	sa->setNumberOfModifications(1);
	TEST_EQUAL (sa->getNumberOfModifications(),1);
	sa->setNumberOfModifications(0);
END_SECTION

START_SECTION(void setTags(const std::vector< String > &tags))
	SuffixArrayTrypticCompressed * satc = new SuffixArrayTrypticCompressed(text,"");
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const vector<String> tags_c (tags);
	satc->setTags(tags);
	vector<String> res = satc->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
END_SECTION

START_SECTION(const std::vector<String>& getTags())
	SuffixArrayTrypticCompressed * satc = new SuffixArrayTrypticCompressed(text,"");
	TEST_EQUAL(satc->getTags().size(),0);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	vector<String> res = satc->getTags();
	TEST_EQUAL(res.at(0),tags.at(0));
	TEST_EQUAL(res.at(1),tags.at(1));
END_SECTION

START_SECTION(void setUseTags(bool use_tags))
	SuffixArrayTrypticCompressed * satc = new SuffixArrayTrypticCompressed(text,"");
	TEST_EQUAL(satc->getUseTags(),false);
	satc->setUseTags(1);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	satc->setUseTags(0);
	TEST_EQUAL(satc->getUseTags(),false);
END_SECTION

START_SECTION(bool getUseTags())
	SuffixArrayTrypticCompressed * satc = new SuffixArrayTrypticCompressed(text,"");
	TEST_EQUAL(satc->getUseTags(),false);
	satc->setUseTags(1);
	TEST_EQUAL(satc->getUseTags(),false);
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const vector<String> tags_c (tags);
	satc->setTags(tags);
	TEST_EQUAL(satc->getUseTags(),true);
	satc->setUseTags(0);
	TEST_EQUAL(satc->getUseTags(),false);
END_SECTION

START_SECTION(bool open(const String &file_name))
	TEST_EXCEPTION (Exception::FileNotFound,sa->open("FileThatNotExists"));
	sa = new SuffixArrayTrypticCompressed(text,"");
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.lcp2"))
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.skip2"))
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.sa2"))		
	sa->save("SuffixArrayTrypticCompressed_test_save");
	SuffixArrayTrypticCompressed * sa2 = new SuffixArrayTrypticCompressed(text,"");
	sa2->open("SuffixArrayTrypticCompressed_test_save");
	TEST_EQUAL(sa->toString(),sa2->toString());
END_SECTION

START_SECTION(bool save(const String &file_name))
	//TEST_EXCEPTION (Exception::UnableToCreateFile,sa->save("/usr/WhereIHaveNoRigths"));
	sa = new SuffixArrayTrypticCompressed(text,"");
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.lcp2"))
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.skip2"))
	NEW_TMP_FILE(String("SuffixArrayTrypticCompressed_test_save.sa2"))
	sa->save("SuffixArrayTrypticCompressed_test_save");
	SuffixArrayTrypticCompressed * sa2 = new SuffixArrayTrypticCompressed(text,"SuffixArrayTrypticCompressed_test_save");
	TEST_EQUAL(sa->toString(),sa2->toString());
END_SECTION

START_SECTION(String toString())
	ptr = new SuffixArrayTrypticCompressed("$","");
	TEST_EQUAL(ptr->toString(),"lcp: 0\nskip: 0");
	String s = sa->toString();
	String sTree = s.substr(0,19);
	String lcp = s.substr(25, 4);
	String skip = s.substr(36,4);
	TEST_EQUAL(sTree,"AA\nAAARAA\nARARP\nARP");
	TEST_EQUAL(lcp,"2120");
	TEST_EQUAL(skip,"1210");
END_SECTION

START_SECTION(void printStatistic())
	NOT_TESTABLE
	//only for internal use
END_SECTION

START_SECTION((void findSpec(std::vector< std::vector< std::pair< std::pair< SignedSize, SignedSize >, double > > > &candidates, const std::vector< double > &spec)))
	double masse[255];
	ResidueDB* rdb = ResidueDB::getInstance();
		
	char aa[] = "ARNDCEQGHILKMFPSTWYV";
		
	for (Size i = 0; i<255;++i)
	{
		masse[i]=0;
	}
	for (Size i = 0; i<strlen(aa);++i)
	{
		const Residue* r = rdb->getResidue(aa[i]);
		masse[(int)aa[i]]=r->getMonoWeight(Residue::Internal);
	}

	sa = new SuffixArrayTrypticCompressed(text, "");
	vector<double> spec;
	//spec.push_back(245.2816);
	spec.push_back(AASequence("AR").getMonoWeight(Residue::Full));
	spec.push_back(AASequence("AAAR").getMonoWeight(Residue::Full));
	//spec.push_back(387.4392);
	vector<double> specc(spec);
	vector<vector<pair<pair<SignedSize, SignedSize>, double> > > res;
	sa->findSpec(res, specc);
	
	TEST_EQUAL(res.size(),specc.size());
	for (Size i = 0; i < res.size(); ++i)
	{
		TEST_EQUAL(res.at(i).size(), 1);
	}
	
	TEST_EQUAL(res.at(0).at(0).first.first, 8)
	TEST_EQUAL(res.at(0).at(0).first.second, 2)
	TEST_EQUAL(res.at(1).at(0).first.first, 1)
	TEST_EQUAL(res.at(1).at(0).first.second, 4)

				
	spec.clear();
	const vector<double> specc2(spec);
	res.clear();
	sa->findSpec(res, specc2);
	TEST_EQUAL(res.size(),0);
	spec.push_back(441.4806);	
	spec.push_back(178.1864);
	const vector<double> specc3 (spec);
	res.clear();
	TEST_EXCEPTION(Exception::InvalidValue, sa->findSpec(res, specc3));
	ifstream i_stream;
	i_stream.open(OPENMS_GET_TEST_DATA_PATH("SuffixArrayTrypticCompressed_test.txt"));
	String txt;
	getline(i_stream,txt);
	
	sa = new SuffixArrayTrypticCompressed(txt,"");
	sa->setNumberOfModifications(0);
  sa->setUseTags(false);
	
	
	vector<double> spec_new;
	for (int i = 500; i < 5000; i += 197)
	{
		spec_new.push_back((double)i);
	}
	const vector<double> specc_new (spec_new);
	res.clear();
	sa->findSpec(res, specc_new);
	//checking for doubled results;
	for (Size i = 0; i < res.size();++i)
	{
		for (Size j = 0;j<res.at(i).size();++j)
		{
			for (Size k = j+1; k < res.at(i).size();++k)
			{
				TEST_EQUAL(res[i][j].first.first==res[i][k].first.first && res[i][j].first.second==res[i][k].first.second, false);
			}
		}
	}

	TOLERANCE_ABSOLUTE(0.55)
	sa->setTolerance(0.5);
		
	// checking if the mass of the found candidates is correct
	// checking if the next character is not a P
	for (Size i = 0; i < res.size();++i)
	{
		for (Size j = 0;j<res.at(i).size();++j)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			double m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length();++k)
			{
				m += masse[(int)seq[k]];
			}
			
			if (txt[res.at(i).at(j).first.first-1]!='$') TEST_NOT_EQUAL(seq[0],'P');
			if (txt[res.at(i).at(j).first.first+res.at(i).at(j).first.second]!='$') TEST_EQUAL(seq[seq.length()-1]=='R'||seq[seq.length()-1]=='K',true)
			
			TEST_REAL_SIMILAR(m,specc_new.at(i));
		}
	}
	// getting all candidates with tags 
	Size number_of_tags=0;
	vector<String> res_with_tags_exp;
	for (Size i = 0; i < res.size();++i)
	{
		for (Size j = 0;j<res.at(i).size();++j)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			bool has_tag = false;
			for (Size k = 2; k < seq.length();++k)
			{
				if (seq.substr(k-2,3)=="AAA"||seq.substr(k-2,3)=="ARA")
				{
					has_tag=true;
					break;
				}
			}
			if (has_tag) 
			{
				++number_of_tags;
				res_with_tags_exp.push_back(seq);
			}
			
		}
	}
	//std::cout<<"number_of_tags_:"<<number_of_tags<<std::endl;
	vector<String> tags;
	tags.push_back("AAA");
	tags.push_back("ARA");
	const vector<String> tags_c (tags);
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
			//if (!has_tag) std::cout <<seq<<std::endl;
			TEST_EQUAL(has_tag, true);
			TEST_EQUAL(res.at(i).at(j).second, 0);
			
			res_with_tags.push_back(seq);
		}
	}
	for (Size i = 0; i < res_with_tags_exp.size();++i)
	{
//		bool was_found = false;
		for (Size j = 0; j < res_with_tags.size();++j)
		{
			if (res_with_tags_exp.at(i)==res_with_tags.at(j))
			{
//				was_found=true;
				break;
			}
		}
		//if (!was_found) //std::cout<<res_with_tags_exp.at(i)<<std::endl;
	}
	//std::cout<<"mod: 1"<<std::endl;
	sa->setNumberOfModifications(1);
	sa->setUseTags(false);
	res.clear();
	sa->findSpec(res, specc_new);
	
	for (Size i = 0; i < res.size();i++)
	{
		for (Size j = 0;j<res.at(i).size();j++)
		{
			String seq = txt.substr(res.at(i).at(j).first.first,res.at(i).at(j).first.second);
			double m = EmpiricalFormula("H2O").getMonoWeight();
			for (Size k = 0; k < seq.length();k++)
			{
				m += masse[(int)seq[k]];
			}
			//if (txt[res.at(i).at(j).first.first+res.at(i).at(j).first.second]=='P')
			//{
				//std::cout<<"hasP:"<<seq<<std::endl;
			//}
			TEST_NOT_EQUAL(txt[res.at(i).at(j).first.first+res.at(i).at(j).first.second],'P');
			TEST_REAL_SIMILAR(m+res.at(i).at(j).second,specc_new.at(i));
			
		}
	}
	// testing if a candidate can belong to several input masses
	spec.clear();
	spec.push_back(441.4806);
	spec.push_back(441.4806);
	const vector<double> specc4 (spec);
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
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
