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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/HANDLERS/MzMLSpectrumDecoder.h>

#define MULTI_LINE_STRING(...) #__VA_ARGS__ 

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(MzMLSpectrumDecoder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzMLSpectrumDecoder* ptr = 0;
MzMLSpectrumDecoder* nullPointer = 0;
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

START_SECTION(( void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
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
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // missing defaultArrayLength -> should give an exception of ParseError
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
}
END_SECTION

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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </NotASpectrum>
  );

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_PRECONDITION_VIOLATED(ptr->domParseSpectrum(testString, cptr))
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // no XML at all here ... 
  ptr = new MzMLSpectrumDecoder();
  std::string testString = MULTI_LINE_STRING(
      Lorem ipsum
  );  

  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // missing 64 bit float tag
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // encode as int instead of float 
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  TEST_EXCEPTION(Exception::ParseError,ptr->domParseSpectrum(testString, cptr))
}
END_SECTION

START_SECTION(([EXTRA] void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // missing m/z array 
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
  );
  OpenMS::Interfaces::SpectrumPtr cptr(new OpenMS::Interfaces::Spectrum);
  ptr->domParseSpectrum(testString, cptr);

  TEST_EQUAL(cptr->getMZArray()->data.size(), 0) // failed since no m/z array is present
  TEST_EQUAL(cptr->getIntensityArray()->data.size(), 0) // failed since no m/z array is present
}
END_SECTION

START_SECTION(( void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // TODO detect multiple occurences of an array
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
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
}
END_SECTION

START_SECTION(( void domParseSpectrum(std::string& in, OpenMS::Interfaces::SpectrumPtr & sptr) ))
{
  // TODO detect multiple occurences of an array
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
            <binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
          </binaryDataArray>

          <binaryDataArray arrayLength="10" encodedLength="108">
            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
            <cvParam cvRef="MS" accession="MS:1000517" name="signal to noise array" />
            <userParam name="name" type="xsd:string" value="binaryDataArray_sn"/>
            <userParam name="name2" type="xsd:string" value="binaryDataArray_sn2"/>
            <binary>AAAAAAAANEAAAAAAAAAyQAAAAAAAADBAAAAAAAAALEAAAAAAAAAoQAAAAAAAACRAAAAAAAAAIEAAAAAAAAAYQAAAAAAAABBAAAAAAAAAAEA=</binary>
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
}
END_SECTION

/// Chromatogram
START_SECTION(( void domParseChromatogram(std::string& in, OpenMS::Interfaces::ChromatogramPtr & sptr) ))
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
            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
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
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

