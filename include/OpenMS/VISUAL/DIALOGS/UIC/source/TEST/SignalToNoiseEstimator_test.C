// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSignalToNoiseEstimator
  : public SignalToNoiseEstimator< >
{
  public:
  TestSignalToNoiseEstimator()
    : SignalToNoiseEstimator< >()
  {
  }

  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf)
  : SignalToNoiseEstimator< >(bpf)
  {
  }

  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
    if (&bpf==this) return *this;

    SignalToNoiseEstimator< >::operator=(bpf);

    return *this;
  }

  protected:

  virtual void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_)
      throw()
  {
    if (scan_first_ == scan_last_)
    {
      std::cout << "bla";
    }
    // do nothing here...
  }

};

START_TEST(SignalToNoiseEstimator, "$Id: SignalToNoiseEstimator_test.C 4855 2009-03-13 01:55:12Z groepl $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = 0;
START_SECTION((SignalToNoiseEstimator()))
	ptr = new TestSignalToNoiseEstimator();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION


START_SECTION((SignalToNoiseEstimator(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  TestSignalToNoiseEstimator sne_copy(sne);
	NOT_TESTABLE
END_SECTION


START_SECTION((SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
	NOT_TESTABLE
END_SECTION


START_SECTION((virtual ~SignalToNoiseEstimator()))
	delete ptr;
END_SECTION


START_SECTION((virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual void init(const Container& c)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual double getSignalToNoise(const PeakIterator& data_point)))
  // hard to do without implementing computeSTN_ properly
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual double getSignalToNoise(const PeakType &data_point)))
  // hard to do without implementing computeSTN_ properly
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


