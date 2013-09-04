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
      <spectrum index="2" id="index=2" spotID="M2" defaultArrayLength="15" dataProcessingRef="dp_sp_2">
        <referenceableParamGroupRef ref="CommonMS1SpectrumParams"/>
        <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="400.39" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="1795.56" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="445.347" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="120054" unitAccession="MS:1000131" unitName="number of counts" unitCvRef="MS"/>
        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="16675500"/>
        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
        <userParam name="sdname" value="spectrumdescription3"/>
        <scanList count="2">
          <cvParam cvRef="MS" accession="MS:1000573" name="median of spectra"/>
          <userParam name="name" value="acquisition_list"/>
          <scan externalSpectrumID="4711">
            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="5.3" unitAccession="UO:0000010" unitName="second" unitCvRef="UO"/>
            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="+ c NSI Full ms [ 400.00-1800.00]"/>
            <cvParam cvRef="MS" accession="MS:1000616" name="preset scan configuration" value="3"/>
            <userParam name="name" value="acquisition1"/>
            <scanWindowList count="1">
              <scanWindow>
                <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="400" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
                <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1800" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              </scanWindow>
            </scanWindowList>
          </scan>
          <scan externalSpectrumID="4712">
            <userParam name="name" value="acquisition2"/>
          </scan>
        </scanList>
        <productList count="1">
          <product>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="18.88"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.0"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="2.0"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <userParam name="iwname" value="isolationwindow3"/>
            </isolationWindow>
          </product>
          <product>
            <isolationWindow>
              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="19.99"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="3.0"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="4.0"  unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
              <userParam name="iwname" value="isolationwindow4"/>
            </isolationWindow>
          </product>
        </productList>
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

