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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/SparseVector.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SparseVector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SparseVector<double>* ptr = nullptr;
SparseVector<double>* nullPointer = nullptr;
START_SECTION(SparseVector())
{
	ptr = new SparseVector<double>();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SparseVector())
{
  delete ptr;
}
END_SECTION

START_SECTION((SparseVector(Value se)))
{
	SparseVector<double> sv(3.0);
	sv.push_back(3.0);
	TEST_EQUAL(sv.size(), 1)
	TEST_EQUAL(sv.nonzero_size(), 0)
}
END_SECTION

SparseVector<double> sv(8,0,3.0);

START_SECTION((SparseVector(size_type size, Value value, Value se=0)))
{
	TEST_EQUAL(sv.size(), 8)
	TEST_EQUAL(sv.nonzero_size(), 8)
}
END_SECTION


SparseVector<double> sv2(sv);

START_SECTION((SparseVector(const SparseVector &source)))
{
	TEST_EQUAL(sv2.size(), 8)
}
END_SECTION

START_SECTION((void resize(size_type newsize)))
{
	sv2.resize(10);
	TEST_EQUAL(sv2.size(),10)
}
END_SECTION

START_SECTION((SparseVector& operator=(const SparseVector &source)))
{
	sv2=sv;
	TEST_EQUAL(sv2.size(), 8)
}
END_SECTION


START_SECTION((bool operator==(const SparseVector &rhs) const ))
{
	SparseVector<double> sv3(sv);
	TEST_EQUAL((sv3==sv), true)
}
END_SECTION

START_SECTION((bool operator<(const SparseVector &rhs) const ))
{
	SparseVector<double> sv3(sv);
	sv3[0]=-1.23;
	TEST_EQUAL((sv3<sv), true)
}
END_SECTION


START_SECTION((void push_back(Value value)))
{
	sv2.push_back(666);
	TEST_EQUAL(sv2.size(), 9)
	TEST_EQUAL(sv2.at(8), 666)
}
END_SECTION

START_SECTION((ValueProxy operator[](size_type pos)))
{
	sv2[8]=3;
	TEST_EQUAL((double)sv2[8],3)
	TEST_EQUAL(sv2.size(), 9)
	TEST_EQUAL(sv2.nonzero_size(), 8)
}
END_SECTION

START_SECTION((const Value operator[](size_type pos) const ))
{
	const SparseVector<double> sv3 = sv2;
	TEST_EQUAL((double)sv3[8],3)
}
END_SECTION


START_SECTION((Value at(size_type pos) const))
{
	TEST_EQUAL(sv2.at(8), 3)
	TEST_EQUAL(sv2.at(0), 0)
}
END_SECTION

START_SECTION((size_type size() const ))
{
	TEST_EQUAL(sv2.size(), 9)
}
END_SECTION

START_SECTION((size_type nonzero_size() const ))
{
	TEST_EQUAL(sv2.nonzero_size(), 8)
}
END_SECTION

START_SECTION((void clear()))
{
	sv2.clear();
	TEST_EQUAL(sv2.size(), 0)
}
END_SECTION

START_SECTION((void erase(SparseVectorIterator it)))
{
	sv.erase(sv.begin()+5);
	TEST_EQUAL(sv.size(),7)

	//real test
	SparseVector<double> sv2;
	sv2.push_back(1.0);
	sv2.push_back(1.1);
	sv2.push_back(1.2);
	sv2.push_back(1.3);
	sv2.push_back(1.4);

	sv2.erase(sv2.begin());
	TEST_EQUAL(sv2.size(),4)
	TEST_EQUAL(sv2.at(0),1.1)
	TEST_EQUAL(sv2.at(1),1.2)
	TEST_EQUAL(sv2.at(2),1.3)
	TEST_EQUAL(sv2.at(3),1.4)

	sv2.erase(sv2.begin()+2);
	TEST_EQUAL(sv2.size(),3)
	TEST_EQUAL(sv2.at(0),1.1)
	TEST_EQUAL(sv2.at(1),1.2)
	TEST_EQUAL(sv2.at(2),1.4)

	sv2.erase(sv2.end()-1);
	TEST_EQUAL(sv2.size(),2)
	TEST_EQUAL(sv2.at(0),1.1)
	TEST_EQUAL(sv2.at(1),1.2)
}
END_SECTION

START_SECTION((void erase(SparseVectorIterator first, SparseVectorIterator last)))
{
	sv[4]=3;
	sv.erase(sv.begin()+5,sv.end());
	TEST_EQUAL(sv.size(),5)

	//real test
	SparseVector<double> sv2;
	sv2.push_back(1.0);
	sv2.push_back(1.1);
	sv2.push_back(1.2);
	sv2.push_back(1.3);
	sv2.push_back(1.4);
	sv2.push_back(1.5);
	sv2.push_back(1.6);
	sv2.push_back(1.7);

	sv2.erase(sv2.begin(),sv2.begin()+2);
	TEST_EQUAL(sv2.size(),6)
	TEST_EQUAL(sv2.at(0),1.2)
	TEST_EQUAL(sv2.at(1),1.3)
	TEST_EQUAL(sv2.at(2),1.4)
	TEST_EQUAL(sv2.at(3),1.5)
	TEST_EQUAL(sv2.at(4),1.6)
	TEST_EQUAL(sv2.at(5),1.7)

	sv2.erase(sv2.begin()+1,sv2.begin()+3);
	TEST_EQUAL(sv2.size(),4)
	TEST_EQUAL(sv2.at(0),1.2)
	TEST_EQUAL(sv2.at(1),1.5)
	TEST_EQUAL(sv2.at(2),1.6)
	TEST_EQUAL(sv2.at(3),1.7)

	sv2.erase(sv2.end()-3,sv2.end());
	TEST_EQUAL(sv2.size(),1)
	TEST_EQUAL(sv2.at(0),1.2)
}
END_SECTION

START_SECTION((SparseVectorIterator getMinElement()))
{
	sv[2]=1;
	sv.erase(sv.begin());
	TEST_EQUAL((double)*(sv.getMinElement()),0)
}
END_SECTION

START_SECTION((iterator begin()))
{
	double i = 0;
	for (SparseVector<double>::iterator vit = sv.begin(); vit != sv.end();++vit)
	{
		i+= (double)*vit;
	}
	TEST_EQUAL(i,4)

	for (SparseVector<double>::iterator vit = sv.end()-1; vit != sv.begin();--vit)
	{
		i-= (double)*vit;
	}
	i-= (double)*(sv.begin());
	TEST_EQUAL(i,0)

	SparseVector<double>::iterator vit = sv.begin();
	vit[1];
	TEST_EQUAL((double)*vit,1)
	vit[2];
	TEST_EQUAL((double)*vit,3)

	vit = sv.begin()+0;
	TEST_EQUAL((double)*vit,0)
	vit = sv.begin()+1;
	TEST_EQUAL((double)*vit,1)
	vit = sv.begin()+2;
	TEST_EQUAL((double)*vit,0)
	vit = sv.begin()+3;
	TEST_EQUAL((double)*vit,3)
	vit -=1;
	TEST_EQUAL((double)*vit,0)
	vit +=1;
	TEST_EQUAL((double)*vit,3)
	vit -=3;
	TEST_EQUAL((double)*vit,0)
	vit +=3;
	TEST_EQUAL((double)*vit,3)
	vit = sv.end()-1;
	TEST_EQUAL((double)*vit,3)
	vit = sv.end()-2;
	TEST_EQUAL((double)*vit,0)
	vit = sv.end()-3;
	TEST_EQUAL((double)*vit,1)
	vit = sv.end()-4;
	TEST_EQUAL((double)*vit,0)

	sv[1] = 3;
	sv[3] = 1;
	vit = sv.begin();
	vit.hop();
	TEST_EQUAL((double)*vit,0)
	vit.hop();
	TEST_EQUAL((double)*vit,1)
	vit.hop();
	//TEST_EQUAL(vit,sv.end())

	TEST_EQUAL(sv.end()-sv.begin(),4)

	TEST_EQUAL(sv.begin()< sv.end(),true)
	TEST_EQUAL(sv.end()> sv.begin(),true)
	TEST_EQUAL(sv.begin()>=sv.begin(),true)
	TEST_EQUAL(sv.end()<= sv.end(),true)
}
END_SECTION

START_SECTION((iterator end()))
{
	NOT_TESTABLE
  // tested above
}
END_SECTION

START_SECTION((const_iterator begin() const ))
{
	double i = 0;
	for (SparseVector<double>::const_iterator cvit = sv.begin(); cvit != sv.end();++cvit)
	{
		i+= (double)*cvit;
	}
	TEST_EQUAL(i,4)

	for (SparseVector<double>::const_iterator cvit = sv.end()-1; cvit != sv.begin();--cvit)
	{
		i-= (double)*cvit;
	}
	i-= (double)*(sv.begin());
	TEST_EQUAL(i,0)

	SparseVector<double>::const_iterator cvit = sv.begin();
	cvit[1];
	TEST_EQUAL((double)*cvit,3)
	cvit[2];
	TEST_EQUAL((double)*cvit,1)

	cvit = sv.begin()+0;
	TEST_EQUAL((double)*cvit,0)
	cvit = sv.begin()+1;
	TEST_EQUAL((double)*cvit,3)
	cvit = sv.begin()+2;
	TEST_EQUAL((double)*cvit,0)
	cvit = sv.begin()+3;
	TEST_EQUAL((double)*cvit,1)
	cvit -=1;
	TEST_EQUAL((double)*cvit,0)
	cvit +=1;
	TEST_EQUAL((double)*cvit,1)
	cvit -=3;
	TEST_EQUAL((double)*cvit,0)
	cvit +=3;
	TEST_EQUAL((double)*cvit,1)
	cvit = sv.end()-1;
	TEST_EQUAL((double)*cvit,1)
	cvit = sv.end()-2;
	TEST_EQUAL((double)*cvit,0)
	cvit = sv.end()-3;
	TEST_EQUAL((double)*cvit,3)
	cvit = sv.end()-4;
	TEST_EQUAL((double)*cvit,0)

	cvit = sv.begin();
	cvit.hop();
	TEST_EQUAL((double)*cvit,0)
	cvit.hop();
	TEST_EQUAL((double)*cvit,1)
	cvit.hop();
	//TEST_EQUAL(cvit,sv.end())

	TEST_EQUAL(sv.end()-sv.begin(),4)

	TEST_EQUAL(sv.begin()< sv.end(),true)
	TEST_EQUAL(sv.end()> sv.begin(),true)
	TEST_EQUAL(sv.begin()>=sv.begin(),true)
	TEST_EQUAL(sv.end()<= sv.end(),true)
}
END_SECTION

START_SECTION((const_iterator end() const ))
{
	NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((reverse_iterator rbegin()))
{
	double i = 0;
	for (SparseVector<double>::reverse_iterator rvit = sv.rbegin(); rvit != sv.rend();++rvit)
	{
		i+= (double)*rvit;
	}
	TEST_EQUAL(i,4)

	for (SparseVector<double>::reverse_iterator rvit = sv.rend()-1; rvit != sv.rbegin();--rvit)
	{
		i-= (double)*rvit;
	}
	i-= (double)*(sv.rbegin());
	TEST_EQUAL(i,0)

	SparseVector<double>::reverse_iterator rvit = sv.rbegin();
	rvit[2];
	TEST_EQUAL((double)*rvit,3)
	rvit[1];
	TEST_EQUAL((double)*rvit,0)

	rvit = sv.rbegin()+0;
	TEST_EQUAL((double)*rvit,1)
	rvit = sv.rbegin()+1;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rbegin()+2;
	TEST_EQUAL((double)*rvit,3)
	rvit = sv.rbegin()+3;
	TEST_EQUAL((double)*rvit,0)
	rvit -=1;
	TEST_EQUAL((double)*rvit,3)
	rvit +=1;
	TEST_EQUAL((double)*rvit,0)
	rvit -=3;
	TEST_EQUAL((double)*rvit,1)
	rvit +=3;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-1;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-2;
	TEST_EQUAL((double)*rvit,3)
	rvit = sv.rend()-3;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-4;
	TEST_EQUAL((double)*rvit,1)


	TEST_EQUAL(sv.rend()-sv.rbegin(),4)

	rvit = sv.rbegin();
	rvit.rhop();
	TEST_EQUAL((double)*rvit,0)
	rvit.rhop();
	TEST_EQUAL((double)*rvit,0)
	rvit.rhop();
	TEST_EQUAL(rvit==sv.rend(),true)
}
END_SECTION

START_SECTION((reverse_iterator rend()))
{
  NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((const_reverse_iterator rbegin() const ))
{
	double i = 0;
	for (SparseVector<double>::const_reverse_iterator rvit = sv.rbegin(); rvit != sv.rend();++rvit)
	{
		i+= (double)*rvit;
	}
	TEST_EQUAL(i,4)

	for (SparseVector<double>::const_reverse_iterator rvit = sv.rend()-1; rvit != sv.rbegin();--rvit)
	{
		i-= (double)*rvit;
	}
	i-= (double)*(sv.rbegin());
	TEST_EQUAL(i,0)

	SparseVector<double>::const_reverse_iterator rvit = sv.rbegin();
	/*
	rvit[2];
	TEST_EQUAL((double)*rvit,3)
	rvit[1];
	TEST_EQUAL((double)*rvit,0)

	rvit = sv.rbegin()+0;
	TEST_EQUAL((double)*rvit,1)
	rvit = sv.rbegin()+1;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rbegin()+2;
	TEST_EQUAL((double)*rvit,3)
	rvit = sv.rbegin()+3;
	TEST_EQUAL((double)*rvit,0)
	rvit -=1;
	TEST_EQUAL((double)*rvit,3)
	rvit +=1;
	TEST_EQUAL((double)*rvit,0)
	rvit -=3;
	TEST_EQUAL((double)*rvit,1)
	rvit +=3;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-1;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-2;
	TEST_EQUAL((double)*rvit,3)
	rvit = sv.rend()-3;
	TEST_EQUAL((double)*rvit,0)
	rvit = sv.rend()-4;
	TEST_EQUAL((double)*rvit,1)
	*/
	rvit = sv.rbegin();
	rvit.rhop();
	TEST_EQUAL((double)*rvit,0)
	rvit.rhop();
	TEST_EQUAL((double)*rvit,0)
	rvit.rhop();
	//TEST_EQUAL(rvit,sv.rend())

	TEST_EQUAL(sv.rend()-sv.rbegin(),4)
}
END_SECTION

START_SECTION((const_reverse_iterator rend() const ))
{
  NOT_TESTABLE
	//tested above
}
END_SECTION

START_SECTION((void print() const))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



