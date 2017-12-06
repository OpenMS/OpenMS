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
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

class IMSAlphabetParserImpl
    : public IMSAlphabetParser<>
{
private:
  ContainerType elements_;
public:
  ContainerType& getElements() override
  {
    return elements_;
  }

  void parse(std::istream& ) override
  {
    // ignore istream, just enter something into the map
    elements_.insert(std::make_pair("A", 71.03711));
    elements_.insert(std::make_pair("R", 156.10111));
  }
};

START_TEST(IMSAlphabetParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// this is an abstract class, that only provides the load method
// it cannot be instanciated so it cannot be tested, therefor we
// test the implementation from above


IMSAlphabetParser<>* ptr = nullptr;
IMSAlphabetParser<>* null_ptr = nullptr;

START_SECTION(IMSAlphabetParser())
{
  ptr = new IMSAlphabetParserImpl();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSAlphabetParser())
{
  delete ptr;
}
END_SECTION

IMSAlphabetParser<> * parser = new IMSAlphabetParserImpl();

START_SECTION((void load(const std::string &fname)))
{
  TEST_EXCEPTION(Exception::IOException ,parser->load(""))

  String filename;
  NEW_TMP_FILE(filename)

  // just create the file
  ofstream of;
  of.open(filename.c_str());
  of << "just text" << std::endl;
  of.close();

  parser->load(filename);

  TEST_EQUAL(parser->getElements().empty(), false)
}
END_SECTION

START_SECTION((virtual ContainerType& getElements()))
{
  TEST_EQUAL(parser->getElements().size(), 2)
  TEST_REAL_SIMILAR(parser->getElements()["A"], 71.03711)
  TEST_REAL_SIMILAR(parser->getElements()["R"], 156.10111)
}
END_SECTION

START_SECTION((virtual void parse(InputSource &is)))
{
  // already tested by load
  NOT_TESTABLE
}
END_SECTION

delete parser;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



