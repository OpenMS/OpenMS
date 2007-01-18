// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/KernelTraits.h>

///////////////////////////

START_TEST(KernelTraits, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

STATUS("KernelTraits struct is a traits class and cannot be instantiated.")
STATUS("Thus we can only employ all its typedef to instantiate something.")
STATUS("But if this test doesn't even COMPILE,")
STATUS("-- well, then there is definitely something going wrong.")

CHECK(KernelTraits)

  KernelTraits::RealType realtype_instance(0);
  KernelTraits::CoordinateType coordinatetype_instance(0);
  KernelTraits::IntensityType intensitytype_instance(0);
  KernelTraits::QualityType qualitytype_instance(0);
  KernelTraits::ProbabilityType probabilitytype_instance(0);
  KernelTraits::PriorityType prioritytype_instance(0);
  KernelTraits::AreaType areatype_instance(0);
  KernelTraits::FullWidthHalfMaxType fullwidthhalfmaxtype_instance(0);
  KernelTraits::WidthType widthtype_instance(0);
  KernelTraits::SignalToNoiseType signaltonoisetype_instance(0);
  KernelTraits::ChargeType chargetype_instance(0);

  realtype_instance++;
  coordinatetype_instance++;
  intensitytype_instance++;
  qualitytype_instance++;
  probabilitytype_instance++;
  prioritytype_instance++;
  areatype_instance++;
  fullwidthhalfmaxtype_instance++;
  widthtype_instance++;
  signaltonoisetype_instance++;
  chargetype_instance++;

RESULT

CHECK(DoubleKernelTraits)

  DoubleKernelTraits::RealType realtype_instance(0);
  DoubleKernelTraits::CoordinateType coordinatetype_instance(0);
  DoubleKernelTraits::IntensityType intensitytype_instance(0);
  DoubleKernelTraits::QualityType qualitytype_instance(0);
  DoubleKernelTraits::ProbabilityType probabilitytype_instance(0);
  DoubleKernelTraits::PriorityType prioritytype_instance(0);
  DoubleKernelTraits::AreaType areatype_instance(0);
  DoubleKernelTraits::FullWidthHalfMaxType fullwidthhalfmaxtype_instance(0);
  DoubleKernelTraits::WidthType widthtype_instance(0);
  DoubleKernelTraits::SignalToNoiseType signaltonoisetype_instance(0);
  DoubleKernelTraits::ChargeType chargetype_instance(0);

  realtype_instance++;
  coordinatetype_instance++;
  intensitytype_instance++;
  qualitytype_instance++;
  probabilitytype_instance++;
  prioritytype_instance++;
  areatype_instance++;
  fullwidthhalfmaxtype_instance++;
  widthtype_instance++;
  signaltonoisetype_instance++;
  chargetype_instance++;

RESULT

CHECK(FloatKernelTraits)

  FloatKernelTraits::RealType realtype_instance(0);
  FloatKernelTraits::CoordinateType coordinatetype_instance(0);
  FloatKernelTraits::IntensityType intensitytype_instance(0);
  FloatKernelTraits::QualityType qualitytype_instance(0);
  FloatKernelTraits::ProbabilityType probabilitytype_instance(0);
  FloatKernelTraits::PriorityType prioritytype_instance(0);
  FloatKernelTraits::AreaType areatype_instance(0);
  FloatKernelTraits::FullWidthHalfMaxType fullwidthhalfmaxtype_instance(0);
  FloatKernelTraits::WidthType widthtype_instance(0);
  FloatKernelTraits::SignalToNoiseType signaltonoisetype_instance(0);
  FloatKernelTraits::ChargeType chargetype_instance(0);

  realtype_instance++;
  coordinatetype_instance++;
  intensitytype_instance++;
  qualitytype_instance++;
  probabilitytype_instance++;
  prioritytype_instance++;
  areatype_instance++;
  fullwidthhalfmaxtype_instance++;
  widthtype_instance++;
  signaltonoisetype_instance++;
  chargetype_instance++;

RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
