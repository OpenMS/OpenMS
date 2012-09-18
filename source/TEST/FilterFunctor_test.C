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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(FilterFunctor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// pure interface, hardly testable

START_SECTION(FilterFunctor())
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor(const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor& operator = (const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(static void registerChildren())
	FilterFunctor* ff = Factory<FilterFunctor>::create("ComplementFilter");
	TEST_EQUAL(ff->getName(), "ComplementFilter")
	ff = Factory<FilterFunctor>::create("IntensityBalanceFilter");
	TEST_EQUAL(ff->getName(), "IntensityBalanceFilter")	
	ff = Factory<FilterFunctor>::create("NeutralLossDiffFilter");
	TEST_EQUAL(ff->getName(), "NeutralLossDiffFilter")
	ff = Factory<FilterFunctor>::create("IsotopeDiffFilter");
	TEST_EQUAL(ff->getName(), "IsotopeDiffFilter")
	ff = Factory<FilterFunctor>::create("TICFilter");
	TEST_EQUAL(ff->getName(), "TICFilter")
END_SECTION

START_SECTION(template<typename SpectrumType> double apply(SpectrumType&))
	NOT_TESTABLE
END_SECTION

START_SECTION(~FilterFunctor())
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
