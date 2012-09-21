// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#ifndef USE_BACKGROUND_CONTROL_H
#define USE_BACKGROUND_CONTROL_H

namespace OpenMS
{

class OPENMS_DLLAPI BackgroundControl
{

	private:

	std::map<double, std::map< double, BackgroundIntensityBin> > intensityBinMap;

	void init();

	public:

	~BackgroundControl();

	BackgroundControl();

	void addPeakMSScan( double , std::list<CentroidPeak>* peakList );

	double getBackgroundLevel( double mz, double tr);

	// find a key in the intensity map:
	std::map<double, std::map< double, BackgroundIntensityBin> >::iterator findTrKey( double );

	// find a key in the m/z map:
	std::map< double, BackgroundIntensityBin>::iterator findMzKey( double mz, std::map< double, BackgroundIntensityBin>* );

	void processIntensityMaps( );

	// overload operators:
	bool operator==(const BackgroundControl&);
	BackgroundControl& operator<=(const BackgroundControl&);
	BackgroundControl& operator>=(const BackgroundControl&);
	BackgroundControl& operator<(const BackgroundControl&);
	BackgroundControl& operator>(const BackgroundControl&);

};

} // namespace

#endif

