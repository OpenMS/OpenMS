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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>

namespace OpenMS
{

	namespace Internal
	{
		template <>
		template <>
		void MzDataHandler <MSExperiment<PickedPeak1D > >::writeDerivedPeakSupplementalData_ < DPeakArray<1, PickedPeak1D > >( std::ostream& os, DPeakArray<1, PickedPeak1D > const & container)
		{
			// default: write data in 32Bit -> fill float array
//			float* tmp = decoder_[0].getFloatBuffer(container.size());
			UInt container_size = container.size();
			// area
//			for (UInt i=0; i<container_size; i++) tmp[i] = container[i].getArea();
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getArea());
			writeBinary_(os,container_size,"supDataArrayBinary","area",1);
			// FWHM
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getFWHM());
			writeBinary_(os,container_size,"supDataArrayBinary","fwhm",2);
			// LeftWidthParameter
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getLeftWidthParameter());
			writeBinary_(os,container_size,"supDataArrayBinary","leftWidth",3);
			// RightWidthParameter
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getRightWidthParameter());
			writeBinary_(os,container_size,"supDataArrayBinary","rightWidth",4);
			// charge
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getCharge());
			writeBinary_(os,container_size,"supDataArrayBinary","charge",5);
			// signal to noise
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getSN());
			writeBinary_(os,container_size,"supDataArrayBinary","signalToNoise",6);
			// rValue
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getRValue());
			writeBinary_(os,container_size,"supDataArrayBinary","rValue",7);
			// peakShape
			data_to_encode_.clear();
			for (UInt i=0; i<container_size; i++) data_to_encode_.push_back (container[i].getPeakShape());
			writeBinary_(os,container_size,"supDataArrayBinary","peakShape",8);
		}

		template <>
		template <>
//		void MzDataHandler <MSExperiment<PickedPeak1D > >::readPeakSupplementalData_ < PickedPeak1D >( std::vector<void*>& data, PickedPeak1D& peak, UInt n)
		void MzDataHandler <MSExperiment<PickedPeak1D > >::readPeakSupplementalData_ < PickedPeak1D >( PickedPeak1D& peak, UInt n)
		{
			enum PickedPeakMembers {AREA = 2, FWHM, LEFT, RIGHT, CHARGE, SN, RVALUE, SHAPE};

			peak.setArea( getDatum_(AREA,n));
			peak.setFWHM( getDatum_(FWHM,n));
			peak.setLeftWidthParameter( getDatum_(LEFT,n));
			peak.setRightWidthParameter( getDatum_(RIGHT,n));
			peak.setCharge(static_cast<Int>(getDatum_(CHARGE,n)));
			peak.setSN( getDatum_(SN,n));
			peak.setRValue( getDatum_(RVALUE,n));
			peak.setPeakShape(PeakShapeType::Enum(int(getDatum_(SHAPE,n))));
/*
			peak.setArea( getDatum_(data,AREA,n));
			peak.setFWHM( getDatum_(data,FWHM,n));
			peak.setLeftWidthParameter( getDatum_(data,LEFT,n));
			peak.setRightWidthParameter( getDatum_(data,RIGHT,n));
			peak.setCharge(static_cast<Int>(getDatum_(data,CHARGE,n)));
			peak.setSN( getDatum_(data,SN,n));
			peak.setRValue( getDatum_(data,RVALUE,n));
			peak.setPeakShape(PeakShapeType::Enum(int(getDatum_(data,SHAPE,n))));
*/
		}
		
	} // namespace Interanal

} // namespace OpenMS


