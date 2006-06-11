// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: MzDataHandler.C,v 1.3 2006/04/24 18:44:10 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>

namespace OpenMS
{

	namespace Internal
	{
		template <>
		template <>
		void MzDataHandler <MSExperiment<DPickedPeak<1,KernelTraits> > >::writeDerivedPeakSupplementalData_ < DPeakArrayNonPolymorphic<1, DPickedPeak<1,KernelTraits> > >( std::ostream& os, DPeakArrayNonPolymorphic<1, DPickedPeak<1,KernelTraits> > const & container)
		{
			// ???? this should work for double as well as float!!!
			float* tmp = decoder_[0].getFloatBuffer(container.size());
			Size container_size = container.size();
			// area
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getArea();
			writeBinary(os,container_size,"supDataArrayBinary","area",1);
			// FWHM
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getFWHM();
			writeBinary(os,container_size,"supDataArrayBinary","fwhm",2);
			// LeftWidthParameter
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getLeftWidthParameter();
			writeBinary(os,container_size,"supDataArrayBinary","leftWidth",3);
			// RightWidthParameter
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getRightWidthParameter();
			writeBinary(os,container_size,"supDataArrayBinary","rightWidth",4);
			// charge
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getCharge();
			writeBinary(os,container_size,"supDataArrayBinary","charge",5);
			// signal to noise
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getSN();
			writeBinary(os,container_size,"supDataArrayBinary","signalToNoise",6);
			// rValue
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getRValue();
			writeBinary(os,container_size,"supDataArrayBinary","rValue",7);
			// peakShape
			for (UnsignedInt i=0; i<container_size; i++) tmp[i] = container[i].getPeakShape();
			writeBinary(os,container_size,"supDataArrayBinary","peakShape",8);
	
		}

	}

} // namespace OpenMS






