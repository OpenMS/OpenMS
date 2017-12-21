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
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/Matrix.h>

using namespace OpenMS;
using namespace std;

START_TEST(TMTTenPlexQuantitationMethod, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TMTTenPlexQuantitationMethod* ptr = nullptr;
TMTTenPlexQuantitationMethod* null_ptr = nullptr;
START_SECTION(TMTTenPlexQuantitationMethod())
{
	ptr = new TMTTenPlexQuantitationMethod();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~TMTTenPlexQuantitationMethod())
{
	delete ptr;
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getName(), "tmt10plex")
}
END_SECTION

START_SECTION((const IsobaricChannelList& getChannelInformation() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = quant_meth.getChannelInformation();

  TEST_EQUAL(channel_list.size(), 10)
  ABORT_IF(channel_list.size() != 10)

  // descriptions are empty by default
  TEST_STRING_EQUAL(channel_list[0].description, "")
  TEST_STRING_EQUAL(channel_list[1].description, "")
  TEST_STRING_EQUAL(channel_list[2].description, "")
  TEST_STRING_EQUAL(channel_list[3].description, "")
  TEST_STRING_EQUAL(channel_list[4].description, "")
  TEST_STRING_EQUAL(channel_list[5].description, "")
  TEST_STRING_EQUAL(channel_list[6].description, "")
  TEST_STRING_EQUAL(channel_list[7].description, "")
  TEST_STRING_EQUAL(channel_list[8].description, "")
  TEST_STRING_EQUAL(channel_list[9].description, "")

  // check masses&co
  TEST_EQUAL(channel_list[0].name, "126")
  TEST_EQUAL(channel_list[0].id, 0)
  TEST_EQUAL(channel_list[0].center, 126.127726)
  TEST_EQUAL(channel_list[0].channel_id_minus_2, -1)
  TEST_EQUAL(channel_list[0].channel_id_minus_1, -1)
  TEST_EQUAL(channel_list[0].channel_id_plus_1, 2)
  TEST_EQUAL(channel_list[0].channel_id_plus_2, 4)

  TEST_EQUAL(channel_list[1].name, "127N")
  TEST_EQUAL(channel_list[1].id, 1)
  TEST_EQUAL(channel_list[1].center, 127.124761)
  TEST_EQUAL(channel_list[1].channel_id_minus_2, -1)
  TEST_EQUAL(channel_list[1].channel_id_minus_1, -1)
  TEST_EQUAL(channel_list[1].channel_id_plus_1, 3)
  TEST_EQUAL(channel_list[1].channel_id_plus_2, 5)

  TEST_EQUAL(channel_list[2].name, "127C")
  TEST_EQUAL(channel_list[2].id, 2)
  TEST_EQUAL(channel_list[2].center, 127.131081)
  TEST_EQUAL(channel_list[2].channel_id_minus_2, -1)
  TEST_EQUAL(channel_list[2].channel_id_minus_1, 0)
  TEST_EQUAL(channel_list[2].channel_id_plus_1, 4)
  TEST_EQUAL(channel_list[2].channel_id_plus_2, 6)

  TEST_EQUAL(channel_list[3].name, "128N")
  TEST_EQUAL(channel_list[3].id, 3)
  TEST_EQUAL(channel_list[3].center, 128.128116)
  TEST_EQUAL(channel_list[3].channel_id_minus_2, -1)
  TEST_EQUAL(channel_list[3].channel_id_minus_1, 1)
  TEST_EQUAL(channel_list[3].channel_id_plus_1, 5)
  TEST_EQUAL(channel_list[3].channel_id_plus_2, 7)

  TEST_EQUAL(channel_list[4].name, "128C")
  TEST_EQUAL(channel_list[4].id, 4)
  TEST_EQUAL(channel_list[4].center, 128.134436)
  TEST_EQUAL(channel_list[4].channel_id_minus_2, 0)
  TEST_EQUAL(channel_list[4].channel_id_minus_1, 2)
  TEST_EQUAL(channel_list[4].channel_id_plus_1, 6)
  TEST_EQUAL(channel_list[4].channel_id_plus_2, 8)

  TEST_EQUAL(channel_list[5].name, "129N")
  TEST_EQUAL(channel_list[5].id, 5)
  TEST_EQUAL(channel_list[5].center, 129.131471)
  TEST_EQUAL(channel_list[5].channel_id_minus_2, 1)
  TEST_EQUAL(channel_list[5].channel_id_minus_1, 3)
  TEST_EQUAL(channel_list[5].channel_id_plus_1, 7)
  TEST_EQUAL(channel_list[5].channel_id_plus_2, 9)

  TEST_EQUAL(channel_list[6].name, "129C")
  TEST_EQUAL(channel_list[6].id, 6)
  TEST_EQUAL(channel_list[6].center, 129.137790)
  TEST_EQUAL(channel_list[6].channel_id_minus_2, 2)
  TEST_EQUAL(channel_list[6].channel_id_minus_1, 4)
  TEST_EQUAL(channel_list[6].channel_id_plus_1, 8)
  TEST_EQUAL(channel_list[6].channel_id_plus_2, -1)

  TEST_EQUAL(channel_list[7].name, "130N")
  TEST_EQUAL(channel_list[7].id, 7)
  TEST_EQUAL(channel_list[7].center, 130.134825)
  TEST_EQUAL(channel_list[7].channel_id_minus_2, 3)
  TEST_EQUAL(channel_list[7].channel_id_minus_1, 5)
  TEST_EQUAL(channel_list[7].channel_id_plus_1, 9)
  TEST_EQUAL(channel_list[7].channel_id_plus_2, -1)

  TEST_EQUAL(channel_list[8].name, "130C")
  TEST_EQUAL(channel_list[8].id, 8)
  TEST_EQUAL(channel_list[8].center, 130.141145)
  TEST_EQUAL(channel_list[8].channel_id_minus_2, 4)
  TEST_EQUAL(channel_list[8].channel_id_minus_1, 6)
  TEST_EQUAL(channel_list[8].channel_id_plus_1, -1)
  TEST_EQUAL(channel_list[8].channel_id_plus_2, -1)

  TEST_EQUAL(channel_list[9].name, "131")
  TEST_EQUAL(channel_list[9].id, 9)
  TEST_EQUAL(channel_list[9].center, 131.138180)
  TEST_EQUAL(channel_list[9].channel_id_minus_2, 5)
  TEST_EQUAL(channel_list[9].channel_id_minus_1, 7)
  TEST_EQUAL(channel_list[9].channel_id_plus_1, -1)
  TEST_EQUAL(channel_list[9].channel_id_plus_2, -1)
}
END_SECTION

START_SECTION((Size getNumberOfChannels() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getNumberOfChannels(), 10)
}
END_SECTION

START_SECTION((virtual Matrix<double> getIsotopeCorrectionMatrix() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;

  // we only check the default matrix here which is an identity matrix
  // for tmt10plex
  Matrix<double> m = quant_meth.getIsotopeCorrectionMatrix();
  TEST_EQUAL(m.rows(), 10)
  TEST_EQUAL(m.cols(), 10)

  ABORT_IF(m.rows() != 10)
  ABORT_IF(m.cols() != 10)

  for(Matrix<double>::SizeType i = 0; i < m.rows(); ++i)
  {
    for(Matrix<double>::SizeType j = 0; j < m.cols(); ++j)
    {
      if (i == j) TEST_REAL_SIMILAR(m(i,j), 1.0)
      else TEST_REAL_SIMILAR(m(i,j), 0.0)
    }
  }
}
END_SECTION

START_SECTION((Size getReferenceChannel() const ))
{
  TMTTenPlexQuantitationMethod quant_meth;
  TEST_EQUAL(quant_meth.getReferenceChannel(), 0)

  Param p;
  p.setValue("reference_channel","128N");
  quant_meth.setParameters(p);

  TEST_EQUAL(quant_meth.getReferenceChannel(), 3)
}
END_SECTION

START_SECTION((TMTTenPlexQuantitationMethod(const TMTTenPlexQuantitationMethod &other)))
{
  TMTTenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "129C");
  qm.setParameters(p);

  TMTTenPlexQuantitationMethod qm2(qm);
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 6)

}
END_SECTION

START_SECTION((TMTTenPlexQuantitationMethod& operator=(const TMTTenPlexQuantitationMethod &rhs)))
{
  TMTTenPlexQuantitationMethod qm;
  Param p = qm.getParameters();
  p.setValue("channel_127N_description", "new_description");
  p.setValue("reference_channel", "130C");
  qm.setParameters(p);

  TMTTenPlexQuantitationMethod qm2 = qm;
  IsobaricQuantitationMethod::IsobaricChannelList channel_list = qm2.getChannelInformation();
  TEST_STRING_EQUAL(channel_list[1].description, "new_description")
  TEST_EQUAL(qm2.getReferenceChannel(), 8)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
