// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIdentificationBase, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(CompNovoIdentificationBase())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIdentificationBase(const CompNovoIdentificationBase &source)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~CompNovoIdentificationBase()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void getIdentifications(std::vector< PeptideIdentification > &ids, const PeakMap &exp)=0))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((CompNovoIdentificationBase& operator=(const CompNovoIdentificationBase &source)))
{
	NOT_TESTABLE
}
END_SECTION


std::set<String>str_set;
str_set.insert("TESTSTRING");
std::set<String>::const_iterator it = str_set.begin();


START_SECTION([CompNovoIdentificationBase::Permut] Permut(const std::set< String >::const_iterator &permut, DoubleReal s))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	TEST_EQUAL(perm.getScore(), 50.0)
	TEST_EQUAL(*perm.getPermut(), "TESTSTRING")
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] Permut(const Permut &rhs))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	CompNovoIdentificationBase::Permut copy(perm);
	TEST_EQUAL(perm.getScore(), copy.getScore())
	TEST_EQUAL(*perm.getPermut(), *copy.getPermut())
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] Permut& operator=(const Permut &rhs))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	CompNovoIdentificationBase::Permut copy(it, 0.0);
	copy=perm;
	TEST_EQUAL(perm.getScore(), copy.getScore())
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] virtual ~Permut())
	CompNovoIdentificationBase::Permut * ptr = new CompNovoIdentificationBase::Permut(it, 50.0);
	delete ptr;
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] void setPermut(const std::set< String >::const_iterator &it))
  std::set<String>str_set;
  str_set.insert("zero");
  std::set<String>::const_iterator it_zero = str_set.begin();
	CompNovoIdentificationBase::Permut perm(it_zero, 50.0);
	perm.setPermut(it);
	TEST_EQUAL(*perm.getPermut(), "TESTSTRING");
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] void setScore(DoubleReal score))
	CompNovoIdentificationBase::Permut perm(it, 50.0);
	perm.setScore(0.0);
	TEST_EQUAL(perm.getScore(), 0.0)
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] DoubleReal getScore() const)
	NOT_TESTABLE //already tested above
END_SECTION

START_SECTION([CompNovoIdentificationBase::Permut] const std::set<String>::const_iterator& getPermut() const)
	NOT_TESTABLE //already tested above
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



