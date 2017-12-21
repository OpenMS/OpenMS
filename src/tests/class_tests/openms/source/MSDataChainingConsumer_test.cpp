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

#include <OpenMS/FORMAT/DATAACCESS/MSDataChainingConsumer.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/NoopMSDataConsumer.h>
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


START_TEST(MSDataChainingConsumer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

MSDataChainingConsumer* chaining_consumer_ptr = nullptr;
MSDataChainingConsumer* chaining_consumer_nullPointer = nullptr;

START_SECTION((MSDataChainingConsumer()))
  chaining_consumer_ptr = new MSDataChainingConsumer();
  TEST_NOT_EQUAL(chaining_consumer_ptr, chaining_consumer_nullPointer)
END_SECTION

START_SECTION((~MSDataChainingConsumer()))
    delete chaining_consumer_ptr;
END_SECTION

START_SECTION(( MSDataChainingConsumer(std::vector<IMSDataConsumer*> consumers) ))
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  chaining_consumer_ptr = new MSDataChainingConsumer(consumer_list);
  TEST_NOT_EQUAL(chaining_consumer_ptr, chaining_consumer_nullPointer)
  delete chaining_consumer_ptr;
END_SECTION

START_SECTION((void consumeSpectrum(SpectrumType & s)))
{
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), true) // nothing happened

  delete chaining_consumer;
}
END_SECTION

START_SECTION(([EXTRA] void consumeSpectrum(SpectrumType & s)))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setSpectraProcessingPtr(FunctionChangeSpectrum);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(transforming_consumer);
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  delete chaining_consumer;

  // note how the transforming consumer still works as deleting the chaining
  // consumer does not take ownership of the consumers
  transforming_consumer->consumeSpectrum(exp.getSpectrum(0) );

  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  delete transforming_consumer;
}
END_SECTION

START_SECTION((void consumeChromatogram(ChromatogramType & c)))
{
  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  chaining_consumer->setExpectedSize(0,1);
  chaining_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), true) // nothing happened

  delete chaining_consumer;
}
END_SECTION

START_SECTION(([EXTRA]void consumeChromatogram(ChromatogramType & c)))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setChromatogramProcessingPtr(FunctionChangeChromatogram);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(transforming_consumer);
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrChromatograms() > 0, true)
  MSChromatogram first_chromatogram = exp.getChromatogram(0);

  chaining_consumer->setExpectedSize(0,1);
  chaining_consumer->consumeChromatogram(exp.getChromatogram(0));

  TEST_EQUAL(first_chromatogram == exp.getChromatogram(0), false) // something happened
  TEST_EQUAL(first_chromatogram.isSorted(), true)
  TEST_EQUAL(exp.getChromatogram(0).isSorted(), false)

  delete chaining_consumer;
}
END_SECTION

START_SECTION((void setExpectedSize(Size, Size)))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((void setExperimentalSettings(const ExperimentalSettings&)))
{
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer();

  chaining_consumer->setExpectedSize(2,0);
  ExperimentalSettings s;
  chaining_consumer->setExperimentalSettings( s );

  TEST_NOT_EQUAL(chaining_consumer, chaining_consumer_nullPointer)
  delete chaining_consumer;
}
END_SECTION

START_SECTION(( void appendConsumer(IMSDataConsumer * consumer) ))
{
  MSDataTransformingConsumer * transforming_consumer = new MSDataTransformingConsumer();
  transforming_consumer->setExpectedSize(2,0);
  transforming_consumer->setSpectraProcessingPtr(FunctionChangeSpectrum);

  std::vector<Interfaces::IMSDataConsumer *> consumer_list;
  consumer_list.push_back(new NoopMSDataConsumer());
  consumer_list.push_back(new NoopMSDataConsumer());
  MSDataChainingConsumer * chaining_consumer = new MSDataChainingConsumer(consumer_list);
  chaining_consumer->appendConsumer(transforming_consumer);

  PeakMap exp;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MzMLFile_1.mzML"), exp);
  TEST_EQUAL(exp.getNrSpectra() > 0, true)
  MSSpectrum first_spectrum = exp.getSpectrum(0);

  chaining_consumer->setExpectedSize(2,0);
  chaining_consumer->consumeSpectrum(exp.getSpectrum(0));

  TEST_EQUAL(first_spectrum == exp.getSpectrum(0), false) // something happened
  TEST_EQUAL(first_spectrum.isSorted(), true)
  TEST_EQUAL(exp.getSpectrum(0).isSorted(), false)

  delete chaining_consumer;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
