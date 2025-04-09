// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

#include <QByteArray>

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

  // Because implementations of zlib and alternatives differ, we just test if 
  // the compressed data requires less space.
  ZlibCompression::compressString(raw_data, compressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_TRUE(compressed_data.size() < raw_data.size())

  ZlibCompression::compressString(raw_data2, compressed_data);
  TEST_EQUAL(raw_data2.size(), 64)
  TEST_TRUE(compressed_data.size() >= raw_data2.size())

  ZlibCompression::compressString(raw_data3, compressed_data);
  TEST_EQUAL(raw_data3.size(), 105)
  TEST_TRUE(compressed_data.size() < raw_data3.size())

  ZlibCompression::compressString(raw_data4, compressed_data);
  TEST_EQUAL(raw_data4.size(), 1052)
  TEST_TRUE(compressed_data.size() < raw_data4.size())
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
  TEST_TRUE(compressed_data.size() < raw_data_q.size())

  ZlibCompression::compressString(raw_data_q2, compressed_data);
  TEST_TRUE(compressed_data.size() >= raw_data_q2.size())

  ZlibCompression::compressString(raw_data_q3, compressed_data);
  TEST_TRUE(compressed_data.size() < raw_data_q3.size())

  ZlibCompression::compressString(raw_data_q4, compressed_data);
  TEST_TRUE(compressed_data.size() < raw_data_q4.size())
}
END_SECTION

START_SECTION((static void uncompressString(const void * compressed_data, size_t nr_bytes, std::string& raw_data)))
{
  std::string compressed_data;
  std::string uncompressed_data;

  ZlibCompression::compressString(raw_data, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data.size(), 58)
  TEST_TRUE(compressed_data.size() < raw_data.size())  
  TEST_EQUAL(uncompressed_data.size(), 58)
  TEST_TRUE(uncompressed_data == raw_data)

  ZlibCompression::compressString(raw_data2, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data2.size(), 64)
  TEST_TRUE(compressed_data.size() >= raw_data2.size())  //  "ABCD..." string is difficult to compress
  TEST_EQUAL(uncompressed_data.size(), 64)
  TEST_TRUE(uncompressed_data == raw_data2)

  ZlibCompression::compressString(raw_data3, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data3.size(), 105)
  TEST_TRUE(compressed_data.size() < raw_data3.size())  
  TEST_EQUAL(uncompressed_data.size(), 105)
  TEST_TRUE(uncompressed_data == raw_data3)

  ZlibCompression::compressString(raw_data4, compressed_data);
  ZlibCompression::uncompressString(&compressed_data[0], compressed_data.size(), uncompressed_data);
  TEST_EQUAL(raw_data4.size(), 1052)
  TEST_TRUE(compressed_data.size() < raw_data4.size())  
  TEST_EQUAL(uncompressed_data.size(), 1052)
  TEST_TRUE(uncompressed_data == raw_data4)
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
  TEST_EQUAL(raw_data_q.size(), 58)
  TEST_TRUE(compressed_data.size() < raw_data_q.size())  
  TEST_EQUAL(uncompressed_data.size(), 58)
  TEST_TRUE(uncompressed_data == raw_data_q)

  ZlibCompression::compressString(raw_data_q2, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q2.size(), 64)
  TEST_TRUE(compressed_data.size() >= raw_data_q2.size())  // difficult to compress...
  TEST_EQUAL(uncompressed_data.size(), 64)
  TEST_TRUE(uncompressed_data == raw_data_q2)

  ZlibCompression::compressString(raw_data_q3, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q3.size(), 105)
  TEST_TRUE(compressed_data.size() < raw_data_q3.size())  
  TEST_EQUAL(uncompressed_data.size(), 105)
  TEST_TRUE(uncompressed_data == raw_data_q3)

  ZlibCompression::compressString(raw_data_q4, compressed_data);
  ZlibCompression::uncompressString(compressed_data, uncompressed_data);
  TEST_EQUAL(raw_data_q4.size(), 1052)
  TEST_TRUE(compressed_data.size() < raw_data_q4.size())  
  TEST_EQUAL(uncompressed_data.size(), 1052)
  TEST_TRUE(uncompressed_data == raw_data_q4)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
