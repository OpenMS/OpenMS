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
// $Maintainer:	David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>


///////////////////////////
#include <OpenMS/ANALYSIS/ID/AScore.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FalseDiscoveryRate, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AScore* ptr = 0;
AScore* nullPointer = 0;
START_SECTION(AScore())
{
	ptr = new AScore();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AScore())
{
	delete ptr;
}
END_SECTION
ptr = new AScore();
START_SECTION((PeptideHit compute(PeptideHit& hit, RichPeakSpectrum& real_spectrum, DoubleReal fmt, Int number_of_phospho_sites)))
{
	
}
END_SECTION
			
START_SECTION((DoubleReal computeCumulativeScore(UInt N,UInt n, DoubleReal p)))
{
	UInt n = 5;
	UInt N = 1;
	DoubleReal p = 0.1;
	DoubleReal score = ptr->computeCumulativeScore(N,n,p);
	TEST_REAL_SIMILAR(score,-1.0);

	n = 1;
	score = ptr->computeCumulativeScore(N,n,p);
	TEST_REAL_SIMILAR(score,0.1);
	N = 3;
	score = ptr->computeCumulativeScore(N,n,p);
	TEST_REAL_SIMILAR(score,0.271);
}
END_SECTION

START_SECTION((void computeHighestPeptides( std::vector< std::vector<DoubleReal> >& peptide_site_scores,std::vector<ProbablePhosphoSites>& sites, vector<vector<Size> >& permutations)))
{
	std::vector< std::vector<DoubleReal> > peptide_site_scores_1;
	std::vector< std::vector<DoubleReal> > peptide_site_scores_2;
	std::vector< std::vector<DoubleReal> > peptide_site_scores_3;
	peptide_site_scores_1.resize(4);
	peptide_site_scores_2.resize(4);
	peptide_site_scores_3.resize(4);
	vector<DoubleReal> temp;
	temp.resize(10);
	for(Size i = 0; i < 10; ++i)
	{
		temp[i] = 0.1;
	}
	peptide_site_scores_1[0] = temp;
	peptide_site_scores_2[3] = temp;
	peptide_site_scores_3[0] = temp;
	temp.clear();
	temp.resize(10);
	for(Size i = 0; i < 10; ++i)
	{
		temp[i] = 0.2;
	}
	peptide_site_scores_1[1] = temp;
	peptide_site_scores_2[0] = temp;
	peptide_site_scores_3[3] = temp;
	temp.clear();
	temp.resize(10);
	for(Size i = 0; i < 10; ++i)
	{
		temp[i] = 0.3;
	}
	peptide_site_scores_1[2] = temp;
	peptide_site_scores_2[1] = temp;
	peptide_site_scores_3[2] = temp;
	temp.clear();
	temp.resize(10);
	for(Size i = 0; i < 10; ++i)
	{
		temp[i] = 0.4;
	}
	peptide_site_scores_1[3] = temp;
	peptide_site_scores_2[2] = temp;
	peptide_site_scores_3[1] = temp;

	
	vector<vector<Size> > permutations;
	vector<Size> per;
	per.push_back(1);
	per.push_back(3);
	per.push_back(5);
	permutations.push_back(per);
	per.clear();
	per.push_back(3);
	per.push_back(5);
	per.push_back(6);	
	permutations.push_back(per);
	per.clear();
	per.push_back(1);
	per.push_back(3);
	per.push_back(6);
	permutations.push_back(per);
	per.clear();
	per.push_back(1);
	per.push_back(5);
	per.push_back(6);
	permutations.push_back(per);

	
	vector<ProbablePhosphoSites> sites;
	ptr->computeHighestPeptides(peptide_site_scores_1,sites,permutations);
	TEST_EQUAL(sites.size(),3)
	TEST_EQUAL(sites[0].seq_1, 3);
	TEST_EQUAL(sites[0].seq_2,1);
	TEST_EQUAL(sites[0].second, 3);
	TEST_EQUAL(sites[0].first,1);	
	TEST_EQUAL(sites[0].peak_depth, 1)
		TEST_EQUAL(sites[1].first,5);
	TEST_EQUAL(sites[1].second,3);
	TEST_EQUAL(sites[1].seq_1, 3);
	TEST_EQUAL(sites[1].seq_2,2);	
	TEST_EQUAL(sites[1].peak_depth, 1)
		TEST_EQUAL(sites[2].first,6);
	TEST_EQUAL(sites[2].second,3);
		TEST_EQUAL(sites[2].seq_1, 3);
	TEST_EQUAL(sites[2].seq_2,0);	
	TEST_EQUAL(sites[2].peak_depth, 1)
	
	ptr->computeHighestPeptides(peptide_site_scores_3,sites,permutations);
	TEST_EQUAL(sites.size(),3)
	TEST_EQUAL(sites[0].seq_1, 1);
	TEST_EQUAL(sites[0].seq_2,3);
	TEST_EQUAL(sites[0].second,1 );
	TEST_EQUAL(sites[0].first,3);	
	TEST_EQUAL(sites[0].peak_depth, 1)
		TEST_EQUAL(sites[1].first,5);
	TEST_EQUAL(sites[1].second,1);
	TEST_EQUAL(sites[1].seq_1, 1);
	TEST_EQUAL(sites[1].seq_2,2);	
	TEST_EQUAL(sites[1].peak_depth, 1)
		TEST_EQUAL(sites[2].first,6);
	TEST_EQUAL(sites[2].second,1);
		TEST_EQUAL(sites[2].seq_1, 1);
	TEST_EQUAL(sites[2].seq_2,0);	
	TEST_EQUAL(sites[2].peak_depth, 1)
	
	ptr->computeHighestPeptides(peptide_site_scores_2,sites,permutations);
	TEST_EQUAL(sites.size(),3)
	TEST_EQUAL(sites[0].seq_1, 2);
	TEST_EQUAL(sites[0].seq_2,1);
	TEST_EQUAL(sites[0].second,5 );
	TEST_EQUAL(sites[0].first,1);	
	TEST_EQUAL(sites[0].peak_depth, 1)
		TEST_EQUAL(sites[1].first,3);
	TEST_EQUAL(sites[1].second,5);
	TEST_EQUAL(sites[1].seq_1, 2);
	TEST_EQUAL(sites[1].seq_2,3);	
	TEST_EQUAL(sites[1].peak_depth, 1)
		TEST_EQUAL(sites[2].first,6);
	TEST_EQUAL(sites[2].second,5);
		TEST_EQUAL(sites[2].seq_1, 2);
	TEST_EQUAL(sites[2].seq_2,0);	
	TEST_EQUAL(sites[2].peak_depth, 1)
	
	peptide_site_scores_1.clear();
	temp.clear();
	temp.resize(10);
	temp[0] = 55;
	temp[1] = 60;
	temp[2]= 75;
	temp[3] = 100;
	temp[4] = 90;
	temp[5] = 120;
	temp[6]  =125;
	temp[7] = 120;
	temp[8] = 100;
	temp[9] = 90;
	peptide_site_scores_1.push_back(temp);
		temp.clear();
	temp.resize(10);
	temp[0] = 40;	
	temp[1] = 50;
	temp[2] = 53;
	temp[3] = 60;
	temp[4] = 50;
	temp[5]= 53;
	temp[6]= 59;
	temp[7] = 53;
	temp[8] = 50;
	temp[9]= 40;
	peptide_site_scores_1.push_back(temp);
	permutations.clear();
	per.clear();
	per.push_back(3);
	permutations.push_back(per);
	per.clear();
	per.push_back(6);
	permutations.push_back(per);
	
	ptr->computeHighestPeptides(peptide_site_scores_1,sites,permutations);
	TEST_EQUAL(sites.size(),1)
	TEST_EQUAL(sites[0].seq_1,0)
	TEST_EQUAL(sites[0].seq_2,1)
	TEST_EQUAL(sites[0].first,3);
	TEST_EQUAL(sites[0].second,6);
	TEST_EQUAL(sites[0].peak_depth, 6)
	
	permutations.clear();
	per.clear();
	per.push_back(3);
	per.push_back(5);
	permutations.push_back(per);
	per.clear();
	per.push_back(5);
	per.push_back(6);
	permutations.push_back(per);
		per.clear();
	per.push_back(3);
	per.push_back(7);
	permutations.push_back(per);
		per.clear();
	per.push_back(3);
	per.push_back(6);
	permutations.push_back(per);
		per.clear();
	per.push_back(5);
	per.push_back(7);
	permutations.push_back(per);
		per.clear();
	per.push_back(6);
	per.push_back(7);
	permutations.push_back(per);
	std::vector< std::vector<DoubleReal> > sec;
	peptide_site_scores_1.push_back(temp);
	peptide_site_scores_1.push_back(temp);
	peptide_site_scores_1.push_back(temp);
	peptide_site_scores_1.push_back(temp);
	ptr->computeHighestPeptides(peptide_site_scores_1,sites,permutations);
	TEST_EQUAL(sites.size(),2)
	TEST_EQUAL(sites[0].seq_1,0)
	TEST_EQUAL(sites[0].seq_2,4)
	TEST_EQUAL(sites[0].first,3);
	TEST_EQUAL(sites[0].second,7);
	TEST_EQUAL(sites[0].peak_depth, 6)	
	TEST_EQUAL(sites[1].seq_1,0)
	TEST_EQUAL(sites[1].seq_2,3)
	TEST_EQUAL(sites[1].first,5);
	TEST_EQUAL(sites[1].second,6);
	TEST_EQUAL(sites[1].peak_depth, 6)	
	
}
END_SECTION

START_SECTION((void compute_site_determining_ions(PeptideHit& hit, Size first,Size second, std::vector<RichPeakSpectrum>& site_determining_ions)))
{
	vector<RichPeakSpectrum> th_spectra;
	RichPeakSpectrum temp1,temp2;
	temp1.setName("VT(Phospho)EQSP");
	temp2.setName("VTEQS(Phospho)P");
	ProbablePhosphoSites candidates;
	candidates.seq_1 = 0;
	candidates.seq_2 = 1;
	candidates.first = 1;
	candidates.second = 4;
	candidates.peak_depth = 1;
	th_spectra.push_back(temp1);
	th_spectra.push_back(temp2);
	vector<RichPeakSpectrum> site_determining_ions;
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),6)
	TEST_EQUAL(site_determining_ions[1].size(),6)
	candidates.first = 4;
	candidates.second = 1;
	candidates.seq_1 = 1;
	candidates.seq_2 = 0;	
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),203.102)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),538.19)
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),201.123)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),540.17)
	
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),6)
	TEST_EQUAL(site_determining_ions[1].size(),6)
	
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),203.102)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),538.19)
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),201.123)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),540.17)

	temp1.setName("T(Phospho)YQYS");
	temp2.setName("TYQYS(Phospho)");
	th_spectra.clear();
	th_spectra.push_back(temp1);
	th_spectra.push_back(temp2);
	candidates.seq_1 = 0;
	candidates.seq_2 = 1;
	candidates.first = 0;
	candidates.second = 4;
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),7)
	TEST_EQUAL(site_determining_ions[1].size(),7)
	
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),106.05)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),636.206)
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),186.016)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),640.201)
	
		candidates.first = 4;
	candidates.second = 0;
	candidates.seq_1 = 1;
	candidates.seq_2 = 0;	
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),7)
	TEST_EQUAL(site_determining_ions[1].size(),7)
	
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),106.05)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),636.206)
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),186.016)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),640.201)	
	
	temp1.setName("TST(Phospho)YQYSYPP");
	temp2.setName("TSTYQYS(Phospho)YPP");
	th_spectra.clear();
	th_spectra.push_back(temp1);
	th_spectra.push_back(temp2);
	candidates.seq_1 = 0;
	candidates.seq_2 = 1;
	candidates.first = 2;
	candidates.second = 6;	
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),9)
	TEST_EQUAL(site_determining_ions[1].size(),9)
	
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),370.101)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),917.403)
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),290.135)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),997.37)

	candidates.seq_1 = 1;
	candidates.seq_2 = 0;
	candidates.first = 6;
	candidates.second = 2;	
	ptr->compute_site_determining_ions(th_spectra,candidates,1,site_determining_ions);
	TEST_EQUAL(site_determining_ions.size(),2)
	TEST_EQUAL(site_determining_ions[0].size(),9)
	TEST_EQUAL(site_determining_ions[1].size(),9)
	
	TEST_REAL_SIMILAR(site_determining_ions[1][0].getMZ(),370.101)
	TEST_REAL_SIMILAR(site_determining_ions[1][site_determining_ions[1].size()-1].getMZ(),917.403)
	TEST_REAL_SIMILAR(site_determining_ions[0][0].getMZ(),290.135)
	TEST_REAL_SIMILAR(site_determining_ions[0][site_determining_ions[0].size()-1].getMZ(),997.37)	
}
END_SECTION

START_SECTION((std::vector<Size> computeTupel_(AASequence& without_phospho)))
AASequence phospho("VTQSPSSP");
vector<Size> tupel(ptr->computeTupel_(phospho));
TEST_EQUAL(4, tupel.size())
TEST_EQUAL(1, tupel[0])
TEST_EQUAL(3,tupel[1])
TEST_EQUAL(5,tupel[2])
TEST_EQUAL(6,tupel[3])
END_SECTION

START_SECTION((std::vector<std::vector<Size> > computePermutations_(std::vector<Size>& tupel,Int number_of_phospho_sites)))
vector<Size> tupel;
tupel.push_back(1);
tupel.push_back(2);
tupel.push_back(3);
tupel.push_back(4);
vector<vector<Size> > permutations;
permutations = ptr->computePermutations_(tupel,1);
TEST_EQUAL(4,permutations.size())
TEST_EQUAL(1,permutations[0][0])
TEST_EQUAL(2,permutations[1][0])
TEST_EQUAL(3,permutations[2][0])
TEST_EQUAL(4,permutations[3][0])

permutations = ptr->computePermutations_(tupel,2);
TEST_EQUAL(6,permutations.size())
TEST_EQUAL(1,permutations[0][0])
TEST_EQUAL(2,permutations[0][1])
TEST_EQUAL(1,permutations[1][0])
TEST_EQUAL(3,permutations[1][1])
TEST_EQUAL(1,permutations[2][0])
TEST_EQUAL(4,permutations[2][1])
TEST_EQUAL(2,permutations[3][0])
TEST_EQUAL(3,permutations[3][1])
TEST_EQUAL(2,permutations[4][0])
TEST_EQUAL(4,permutations[4][1])
TEST_EQUAL(3,permutations[5][0])
TEST_EQUAL(4,permutations[5][1])

permutations = ptr->computePermutations_(tupel,3);

TEST_EQUAL(4,permutations.size())
TEST_EQUAL(1,permutations[0][0])
TEST_EQUAL(2,permutations[0][1])
TEST_EQUAL(3,permutations[0][2])
TEST_EQUAL(1,permutations[1][0])
TEST_EQUAL(2,permutations[1][1])
TEST_EQUAL(4,permutations[1][2])
TEST_EQUAL(1,permutations[2][0])
TEST_EQUAL(3,permutations[2][1])
TEST_EQUAL(4,permutations[2][2])
TEST_EQUAL(2,permutations[3][0])
TEST_EQUAL(3,permutations[3][1])
TEST_EQUAL(4,permutations[3][2])


permutations = ptr->computePermutations_(tupel,4);
TEST_EQUAL(1,permutations.size())
TEST_EQUAL(1,permutations[0][0])
TEST_EQUAL(2,permutations[0][1])
TEST_EQUAL(3,permutations[0][2])
TEST_EQUAL(4,permutations[0][3])

END_SECTION
delete ptr;
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
