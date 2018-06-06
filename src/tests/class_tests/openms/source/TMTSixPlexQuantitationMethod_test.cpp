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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTSixPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TMTSixPlexQuantitationMethod* ptr = nullptr;
TMTSixPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(TMTSixPlexQuantitationMethod())
{
	ptr = new TMTSixPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTSixPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getName(), "tmt6plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();
  
  TEST_EQUAL(channel_list.size(), 6)
  ABORT_IF(channel_list.size() != 6)
  
  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
  TEST_STRING_EQUAL(channel_list[4].description, "")
  TEST_STRING_EQUAL(channel_list[5].description, "")
    
  // check masses&co
  TEST_EQUAL(channel_list[0].name, 126)
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 126.127725)
        
  TEST_EQUAL(channel_list[1].name, 127)
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 127.124760)

  TEST_EQUAL(channel_list[2].name, 128)
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 128.134433)

  TEST_EQUAL(channel_list[3].name, 129)
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 129.131468)

  TEST_EQUAL(channel_list[4].name, 130)
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 130.141141)

  TEST_EQUAL(channel_list[5].name, 131)
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 131.138176)

}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 6)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  
  // we only check the default matrix here which is an identity matrix 
  // for tmt6plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 6)
  TEST_EQUAL(m.cols(), 6)
    
  ABORT_IF(m.rows() != 6)
  ABORT_IF(m.cols() != 6)
  
  for(Matrix<double>::SizeType i = 0; i < m.rows(); ++i)
  {
    for(Matrix<double>::SizeType j = 0; j < m.cols(); ++j)
    {
      if (i == j)
        TEST_REAL_SIMILAR(m(i,j), 1.0)
      else
        TEST_REAL_SIMILAR(m(i,j), 0.0)
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTSixPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)
    
  Param p;
  p.setValue("reference_channel",128);
  quant_meth.setParameters(p);
  
  TEST_EQUAL(quant_meth.getReferenceChannel(), 2)
}
END_SECTION

START_SECTION((TMTSixPlexQuantitationMethod(const TMTSixPlexQuantitationMethod &other)))
{
  TMTSixPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127_description", "new_description");
  p.setValue("reference_channel", 129);
  qm.setParameters(p);
  
  TMTSixPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)

}
END_SECTION

START_SECTION((TMTSixPlexQuantitationMethod& operator=(const TMTSixPlexQuantitationMethod &rhs)))
{
  TMTSixPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127_description", "new_description");
  p.setValue("reference_channel", 129);
  qm.setParameters(p);
  
  TMTSixPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
