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

///////////////////////////
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SuffixArray, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((SuffixArray()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((SuffixArray(const String &st, const String &filename)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((SuffixArray(const SuffixArray &sa)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual ~SuffixArray()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual String toString()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void findSpec(std::vector< std::vector< std::pair< std::pair< SignedSize, SignedSize >, DoubleReal > > > &candidates, const std::vector< DoubleReal > &spec)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool save(const String &filename)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool open(const String &filename)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setTolerance(DoubleReal t)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual DoubleReal getTolerance() const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool isDigestingEnd(const char aa1, const char aa2) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setTags(const std::vector< String > &tags)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual const std::vector<String>& getTags()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setUseTags(bool use_tags)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual bool getUseTags()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void setNumberOfModifications(Size number_of_mods)=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual Size getNumberOfModifications()=0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void printStatistic()=0))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



