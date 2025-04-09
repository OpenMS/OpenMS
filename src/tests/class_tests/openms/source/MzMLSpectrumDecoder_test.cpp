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

#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>

#define MULTI_LINE_STRING(...) #__VA_ARGS__ 

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(MzMLSpectrumDecoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLSpectrumDecoder* ptr = nullptr;
MzMLSpectrumDecoder* nullPointer = nullptr;
START_SECTION((MzMLSpectrumDecoder() ))
{
  ptr = new MzMLSpectrumDecoder();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~MzMLSpectrumDecoder()))
{
  delete ptr;
}
END_SECTION

// Working example of parsing a spectrum
START_SECTION(( void domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
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

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 15)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 15)

  TEST_REAL_SIMILAR(cptr->getMZArray()->data[7], 7)
  TEST_REAL_SIMILAR(cptr->getIntensityArray()->data[7], 8)
  delete ptr;
}
END_SECTION

// Working example of parsing a spectrum with some extra CV terms in there (should still work)
START_SECTION(([EXTRA] void domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 15)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 15)

  TEST_REAL_SIMILAR(cptr->getMZArray()->data[7], 7)
  TEST_REAL_SIMILAR(cptr->getIntensityArray()->data[7], 8)
  delete ptr;
}
END_SECTION

// missing defaultArrayLength -> should give an exception of ParseError
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
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

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// root tag is neither spectrum or chromatogram -> precondition violation
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // root tag is neither spectrum or chromatogram
  //
  // this does not generate a runtime error but rather a precondition violation
  // -> it should allow a developer to easily spot a problem with the code if
  // some other XML tag is used.
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <NotASpectrum index="2" id="index=2">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </NotASpectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_PRECONDITION_VIOLATED(ptr->domParseSpectrum(testString, cptr))

  // TEST_PRECONDITION_VIOLATED should already be sufficient to not trigger the "no subtests performed in"
  TEST_EQUAL(cptr->getMZArray()->data.size(), 0)
  delete ptr;
}
END_SECTION

// no XML at all here ...  -> Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      Lorem ipsum
  );  

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// Working example without intensity -> simply an empty spectrum
START_SECTION(( void domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 0)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 0)
  delete ptr;
}
END_SECTION

// missing 64 bit float tag -> should throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
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
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// This is a valid XML structure, but simply empty <binary></binary> -> empty spectra output
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary></binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary></binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 0)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 0)
  delete ptr;
}
END_SECTION

// Invalid XML (unclosed brackets) -> should throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <bina
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// Invalid XML (unclosed brackets) -> should throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <cvParam cvRef="MS" accession="MS:100057"
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// Invalid mzML (too much content inside <binary>) -> should throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>
              <whoPutMeHere>
                some crazy person, obviously
              </whoPutMeHere>
            </binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// Invalid mzML (missing <binary> tag)-> should throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// Invalid content of <binary> -> empty spectrum
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="0" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>
              whoPutMeHere: some crazy person, obviously! What if I contain invalid characters like these &- 
            </binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ConversionError, ptr->domParseSpectrum(testString, cptr) );
  delete ptr;
}
END_SECTION

// encode as int instead of float -> throw Exception
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000519" name="32-bit int" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
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
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
  delete ptr;
}
END_SECTION

// missing m/z array -> no Exception but simply empty data
START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
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
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 0) // failed since no m/z array is present
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 0) // failed since no m/z array is present
  delete ptr;
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // missing: detect semantically invalid XML structures
  // for example: multiple occurrences of an array
  // (fix in MzMLHandlerHelper::computeDataProperties_)
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 15)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 15)

  TEST_REAL_SIMILAR(cptr->getMZArray()->data[7], 7)
  TEST_REAL_SIMILAR(cptr->getIntensityArray()->data[7], 8)
  delete ptr;
}
END_SECTION

// Working example of parsing a chromatogram
START_SECTION(( void domParseChromatogram(const std::string& in, OpenMS::Interfaces::ChromatogramPtr & cptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING( 
      <chromatogram index="1" id="sic native" defaultArrayLength="10" >
        <cvParam cvRef="MS" accession="MS:1000235" name="total ion current chromatogram" value=""/>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="108" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000595" name="time array" unitAccession="UO:0000010" unitName="second" unitCvRef="UO"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkA=</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="108" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary>AAAAAAAAJEAAAAAAAAAiQAAAAAAAACBAAAAAAAAAHEAAAAAAAAAYQAAAAAAAABRAAAAAAAAAEEAAAAAAAAAIQAAAAAAAAABAAAAAAAAA8D8=</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </chromatogram>);

  OpenMS::Interfaces::ChromatogramPtr cptr(new OpenMS::Interfaces::Chromatogram);
  ptr->domParseChromatogram(testString, cptr);

  TEST_EQUAL(cptr->getTimeArray()->data.size(), 10)
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 10)

  TEST_REAL_SIMILAR(cptr->getTimeArray()->data[5], 5)
  TEST_REAL_SIMILAR(cptr->getIntensityArray()->data[5], 5)
  delete ptr;
}
END_SECTION

START_SECTION(( void domParseSpectrum(const std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      <spectrum index="2" id="index=2" defaultArrayLength="15">
        <binaryDataArrayList count="3">
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000786" name="non-standard data array" value="Ion Mobility" />
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  MSSpectrum s;
  ptr->domParseSpectrum(testString, s);

  TEST_EQUAL(s.size(), 15)
  TEST_EQUAL(s.getFloatDataArrays().size(), 1)

  TEST_REAL_SIMILAR(s[7].getMZ(), 7)
  TEST_REAL_SIMILAR(s[7].getIntensity(), 8)
  TEST_REAL_SIMILAR(s.getFloatDataArrays()[0][7], 8)
  TEST_EQUAL(s.getFloatDataArrays()[0].getName(), "Ion Mobility")
  delete ptr;
}
END_SECTION

START_SECTION(( void domParseChromatogram(const std::string& in, OpenMS::Interfaces::ChromatogramPtr & cptr) ))
{
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING( 
      <chromatogram index="1" id="sic native" defaultArrayLength="10" >
        <cvParam cvRef="MS" accession="MS:1000235" name="total ion current chromatogram" value=""/>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="108" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000595" name="time array" unitAccession="UO:0000010" unitName="second" unitCvRef="UO"/>
            <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkA=</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="108" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
            <binary>AAAAAAAAJEAAAAAAAAAiQAAAAAAAACBAAAAAAAAAHEAAAAAAAAAYQAAAAAAAABRAAAAAAAAAEEAAAAAAAAAIQAAAAAAAAABAAAAAAAAA8D8=</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="160" >
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
            <cvParam cvRef="MS" accession="MS:1000786" name="non-standard data array" value="Ion Mobility" />
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </chromatogram>);

  MSChromatogram s;
  ptr->domParseChromatogram(testString, s);

  TEST_EQUAL(s.size(), 10)
  TEST_EQUAL(s.getFloatDataArrays().size(), 1)

  TEST_REAL_SIMILAR(s[5].getRT(), 5)
  TEST_REAL_SIMILAR(s[5].getIntensity(), 5)
  TEST_REAL_SIMILAR(s.getFloatDataArrays()[0][7], 8)
  TEST_EQUAL(s.getFloatDataArrays()[0].getName(), "Ion Mobility")
  delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

