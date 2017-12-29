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
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

void FunctionChangeSpectrum (OpenMS::MSSpectrum & s)
{
  s.sortByIntensity();
}

void FunctionChangeChromatogram (OpenMS::MSChromatogram & c)
{
  c.sortByIntensity();
}


START_TEST(MSDataTransformingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataTransformingConsumer* transforming_consumer_ptr = nullptr;
MSDataTransformingConsumer* transforming_consumer_nullPointer = nullptr;

START_SECTION((MSDataTransformingConsumer()))
  transforming_consumer_ptr = new MSDataTransformingConsumer();
  TEST_NOT_EQUAL(transforming_consumer_ptr, transforming_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataTransformingConsumer()))
    delete transforming_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), true) // nothing happened

  delete transforming_consumer;
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  transforming_consumer->setExpectedSize(0,1);
  transforming_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), true) // nothing happened

  delete transforming_consumer;
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

  transforming_consumer->setExpectedSize(2,0);
  ExperimentalSettings s;
  transforming_consumer->setExperimentalSettings( s );

  TEST_NOT_EQUAL(transforming_consumer, transforming_consumer_nullPointer)
  delete transforming_consumer;
}
END_SECTION

START_SECTION((virtual void setSpectraProcessingPtr( void (*sproptr)(SpectrumType&) )))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  exp.getSpectrum(0).sortByPosition();
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setSpectraProcessingPtr(FunctionChangeSpectrum);
  transforming_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  delete transforming_consumer;
}
END_SECTION

START_SECTION((virtual void setChromatogramProcessingPtr( void (*cproptr)(ChromatogramType&) )))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  exp.getChromatogram(0).sortByPosition();
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  transforming_consumer->setExpectedSize(0,1);
  transforming_consumer->setChromatogramProcessingPtr(FunctionChangeChromatogram);
  transforming_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), false) // something happened
  TEST_EQUAL(first_chromatogram.isSorted(), true)
  TEST_EQUAL(exp.getChromatogram(0).isSorted(), false)

  delete transforming_consumer;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
