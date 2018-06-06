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
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/PIP/LocalLinearMap.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(LocalLinearMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LocalLinearMap* ptr = nullptr;
LocalLinearMap* nullPointer = nullptr;
LocalLinearMap llm;


START_SECTION(LocalLinearMap())
	ptr = new LocalLinearMap();
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getLLMParam().xdim, llm.getLLMParam().xdim)
	TEST_EQUAL(ptr->getLLMParam().ydim, llm.getLLMParam().ydim)
	TEST_EQUAL(ptr->getLLMParam().radius, llm.getLLMParam().radius)
END_SECTION

START_SECTION((virtual ~LocalLinearMap()))
	delete ptr;
END_SECTION

START_SECTION((const LLMParam& getLLMParam() const))
	TEST_EQUAL(llm.getLLMParam().xdim, 1)
	TEST_EQUAL(llm.getLLMParam().ydim, 2)
	TEST_EQUAL(llm.getLLMParam().radius, 0.4)
END_SECTION

START_SECTION((const Matrix<double>& getCodebooks() const))
	TEST_EQUAL(llm.getCodebooks().rows(), 2)
	TEST_EQUAL(llm.getCodebooks().cols(), 18)
	//-0.06281751 0.9460272
	//0.03852812 -0.4956029
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(0,0), -0.06281751)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(0,1), 0.9460272)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(1,0), 0.03852812)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(1,1), -0.4956029)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(0,17), 0.3478902)
	TEST_REAL_SIMILAR(llm.getCodebooks().getValue(1,17), -0.1460901)
END_SECTION

START_SECTION((const Matrix<double>& getMatrixA() const))
	TEST_EQUAL(llm.getMatrixA().rows(), 2)
	TEST_EQUAL(llm.getMatrixA().cols(), 18)
	//-0.005066359 -0.0251465
	//-0.221425369 -0.2565968
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(0,0), -0.005066359)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(0,1), -0.0251465)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(1,0), -0.221425369)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(1,1), -0.2565968)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(0,17), -0.3692879)
	TEST_REAL_SIMILAR(llm.getMatrixA().getValue(1,17), 0.3665653)
END_SECTION

START_SECTION((const vector<double>& getVectorWout() const))
	TEST_EQUAL(llm.getVectorWout().size(), 2)
	//3.746677 
	//3.395571
	TEST_REAL_SIMILAR(llm.getVectorWout()[0], 3.746677)
	TEST_REAL_SIMILAR(llm.getVectorWout()[1], 3.395571)
END_SECTION

START_SECTION((const Matrix<UInt>& getCord() const))
	TEST_EQUAL(llm.getCord().rows(), 2)
	TEST_EQUAL(llm.getCord().cols(), 2)
	TEST_EQUAL(llm.getCord().getValue(0,0), 0)
	TEST_EQUAL(llm.getCord().getValue(0,1), 0)
	TEST_EQUAL(llm.getCord().getValue(1,0), 0)
	TEST_EQUAL(llm.getCord().getValue(1,1), 1)
END_SECTION

START_SECTION((std::vector<double> neigh(const Matrix< UInt > &cord, Size win, double radius)))
	vector<double> nei1 = llm.neigh(llm.getCord(), 0, llm.getLLMParam().radius);
	TEST_EQUAL(nei1[0], 1)
	TEST_REAL_SIMILAR(nei1[1], 0.04393693) 	
END_SECTION

START_SECTION((void normalizeVector(std::vector< double > &aaIndexVariables)))
	NOT_TESTABLE
END_SECTION


END_TEST
