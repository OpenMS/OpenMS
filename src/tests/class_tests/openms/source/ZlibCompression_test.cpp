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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ZlibCompression.h>
///////////////////////////

#define MULTI_LINE_STRING(...) #__VA_ARGS__ 


START_TEST(ZlibCompression, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

ZlibCompression* zlib_ptr = nullptr;

START_SECTION((ZlibCompression()))
  zlib_ptr = new ZlibCompression();
END_SECTION

START_SECTION((~ZlibCompression()))
  delete zlib_ptr;
END_SECTION

char const * toCompress = "AAAAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
char const * toCompress2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
char const * toCompress3 = "Freude, schoner Gotterfunken, Tochter aus Elysium, Wir betreten feuertrunken, Himmlische, dein Heiligtum!";

std::string raw_data(toCompress);
std::string raw_data2(toCompress2);
std::string raw_data3(toCompress3);
std::string raw_data4 = MULTI_LINE_STRING(
    <spectrum index="2" id="index=2" defaultArrayLength="15">
      <binaryDataArrayList count="2">
        <binaryDataArray encodedLength="160" >
          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
          <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
          <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
        </binaryDataArray>
        <binaryDataArray encodedLength="160" >
          <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
          <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
          <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
          <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
        </binaryDataArray>
      </binaryDataArrayList>
    </spectrum>
);

START_SECTION((static void compressString(std::string& raw_data, std::string& compressed_data)))
{
  std::string compressed_data;

  ZlibCompression::compressString(raw_data, compressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_EQUAL(compressed_data.size(), 14)

  ZlibCompression::compressString(raw_data2, compressed_data);
  TEST_EQUAL(raw_data2.size(), 64)
  TEST_EQUAL(compressed_data.size(), 72)

  ZlibCompression::compressString(raw_data3, compressed_data);
  TEST_EQUAL(raw_data3.size(), 105)
  TEST_EQUAL(compressed_data.size(), 97)

  ZlibCompression::compressString(raw_data4, compressed_data);
  TEST_EQUAL(raw_data4.size(), 1052)
  TEST_EQUAL(compressed_data.size(), 335)
}
END_SECTION

START_SECTION((static void compressString(const QByteArray& raw_data, QByteArray& compressed_data)))
{
  QByteArray raw_data_q = QByteArray::fromRawData(&raw_data[0], raw_data.size());
  QByteArray raw_data_q2 = QByteArray::fromRawData(&raw_data2[0], raw_data2.size());
  QByteArray raw_data_q3 = QByteArray::fromRawData(&raw_data3[0], raw_data3.size());
  QByteArray raw_data_q4 = QByteArray::fromRawData(&raw_data4[0], raw_data4.size());
  QByteArray compressed_data;

  ZlibCompression::compressString(raw_data_q, compressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_EQUAL(compressed_data.size(), 14)

  ZlibCompression::compressString(raw_data_q2, compressed_data);
  TEST_EQUAL(raw_data_q2.size(), 64)
  TEST_EQUAL(compressed_data.size(), 72)

  ZlibCompression::compressString(raw_data_q3, compressed_data);
  TEST_EQUAL(raw_data_q3.size(), 105)
  TEST_EQUAL(compressed_data.size(), 97)

  ZlibCompression::compressString(raw_data_q4, compressed_data);
  TEST_EQUAL(raw_data_q4.size(), 1052)
  TEST_EQUAL(compressed_data.size(), 335)
}
END_SECTION

START_SECTION((static void uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data)))
{
  std::string compressed_data;
  std::string uncompressed_data;

  ZlibCompression::compressString(raw_data, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_EQUAL(compressed_data.size(), 14)
  TEST_EQUAL(uncompressed_data.size(), 58)
  TEST_EQUAL(uncompressed_data == raw_data, true)

  ZlibCompression::compressString(raw_data2, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data2.size(), 64)
  TEST_EQUAL(compressed_data.size(), 72)
  TEST_EQUAL(uncompressed_data.size(), 64)
  TEST_EQUAL(uncompressed_data == raw_data2, true)

  ZlibCompression::compressString(raw_data3, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data3.size(), 105)
  TEST_EQUAL(compressed_data.size(), 97)
  TEST_EQUAL(uncompressed_data.size(), 105)
  TEST_EQUAL(uncompressed_data == raw_data3, true)

  ZlibCompression::compressString(raw_data4, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data4.size(), 1052)
  TEST_EQUAL(compressed_data.size(), 335)
  TEST_EQUAL(uncompressed_data.size(), 1052)
  TEST_EQUAL(uncompressed_data == raw_data4, true)
}
END_SECTION
  
START_SECTION((static void uncompressString(const QByteArray& compressed_data, QByteArray& raw_data)))
{
  QByteArray raw_data_q = QByteArray::fromRawData(&raw_data[0], raw_data.size());
  QByteArray raw_data_q2 = QByteArray::fromRawData(&raw_data2[0], raw_data2.size());
  QByteArray raw_data_q3 = QByteArray::fromRawData(&raw_data3[0], raw_data3.size());
  QByteArray raw_data_q4 = QByteArray::fromRawData(&raw_data4[0], raw_data4.size());

  QByteArray compressed_data;
  QByteArray uncompressed_data;

  ZlibCompression::compressString(raw_data_q, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_EQUAL(compressed_data.size(), 14)
  TEST_EQUAL(uncompressed_data.size(), 58)
  TEST_EQUAL(uncompressed_data == raw_data_q, true)

  ZlibCompression::compressString(raw_data_q2, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q2.size(), 64)
  TEST_EQUAL(compressed_data.size(), 72)
  TEST_EQUAL(uncompressed_data.size(), 64)
  TEST_EQUAL(uncompressed_data == raw_data_q2, true)

  ZlibCompression::compressString(raw_data_q3, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q3.size(), 105)
  TEST_EQUAL(compressed_data.size(), 97)
  TEST_EQUAL(uncompressed_data.size(), 105)
  TEST_EQUAL(uncompressed_data == raw_data_q3, true)

  ZlibCompression::compressString(raw_data_q4, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q4.size(), 1052)
  TEST_EQUAL(compressed_data.size(), 335)
  TEST_EQUAL(uncompressed_data.size(), 1052)
  TEST_EQUAL(uncompressed_data == raw_data_q4, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
