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
// $Maintainer: Eva Lange $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>


#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////

START_TEST(LinearResampler, "$Id: LinearResampler_test.C 6084 2009-10-06 00:34:12Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

LinearResampler* lr_ptr = 0;
START_SECTION((LinearResampler()))
  lr_ptr = new LinearResampler;
  TEST_NOT_EQUAL(lr_ptr,0);
END_SECTION

START_SECTION((~LinearResampler()))
  delete lr_ptr;
END_SECTION

START_SECTION((template<typename PeakType> void raster(MSSpectrum<PeakType>& spectrum)))
  MSSpectrum< Peak1D > spec;
  spec.resize(5);
  spec[0].setMZ(0);
  spec[0].setIntensity(3.0f);
  spec[1].setMZ(0.5);
  spec[1].setIntensity(6.0f);
  spec[2].setMZ(1.);
  spec[2].setIntensity(8.0f);
  spec[3].setMZ(1.6);
  spec[3].setIntensity(2.0f);
  spec[4].setMZ(1.8);
  spec[4].setIntensity(1.0f);

	LinearResampler lr;
	Param param;
	param.setValue("spacing",0.5);
  lr.setParameters(param);
  lr.raster(spec);

  DoubleReal sum = 0.0;
	for (Size i=0; i<spec.size(); ++i)
  {
    sum += spec[i].getIntensity();
  }
  TEST_REAL_SIMILAR(sum, 20);
END_SECTION

START_SECTION(( template <typename PeakType > void rasterExperiment(MSExperiment<PeakType>& exp)))
  MSSpectrum< RichPeak1D > spec;
  spec.resize(5);
  spec[0].setMZ(0);
  spec[0].setIntensity(3.0f);
  spec[1].setMZ(0.5);
  spec[1].setIntensity(6.0f);
  spec[2].setMZ(1.);
  spec[2].setIntensity(8.0f);
  spec[3].setMZ(1.6);
  spec[3].setIntensity(2.0f);
  spec[4].setMZ(1.8);
  spec[4].setIntensity(1.0f);

  MSExperiment< RichPeak1D > exp;
  exp.push_back(spec);
  exp.push_back(spec);

  LinearResampler lr;
  Param param;
	param.setValue("spacing",0.5);
  lr.setParameters(param);
  lr.rasterExperiment(exp);


	for (Size s=0; s<exp.size(); ++s)
  {
	  DoubleReal sum = 0.0;
		for (Size i=0; i<exp[s].size(); ++i)
	  {
	    sum += exp[s][i].getIntensity();
	  }
	  TEST_REAL_SIMILAR(sum, 20);
	}

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
