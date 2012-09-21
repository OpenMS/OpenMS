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

#ifndef USE_BACKGROUND_INTENSITY_BIN_H
#define USE_BACKGROUND_INTENSITY_BIN_H

namespace OpenMS
{

class OPENMS_DLLAPI BackgroundIntensityBin
{

	private:

	// do not allow default constructor
	BackgroundIntensityBin()
	{
	};

	// mz and tr coordinates of the bin:
	double mzCoord_;
	double trCoord_;
	double zCoord_;

	std::vector<double> intensityMap_;
	std::map<double, double> intensityHist_;

	double mean_;
	void computeIntensityHist();

	public:

	virtual ~BackgroundIntensityBin();

	// copy constructor
	BackgroundIntensityBin(const BackgroundIntensityBin& );

	// assignment operator
	BackgroundIntensityBin& operator = (const BackgroundIntensityBin& );

	BackgroundIntensityBin(double, double);

	/*
	 * @brief Check if a peak belongs to this intenisty bin
	 */
	bool checkBelonging(MSPeak*);

	/*
	 * @brief Add intensity to BackgroundIntensityBin
	 */
	void addIntensity( double );

	/*
	 * @brief add peak to BackgroundIntensityBin
	 */
	void addMSPeak( MSPeak* );
	/*
	 * @brief Process collected intensities in the map
	 */
	void processIntensities( );

	std::vector<double>* getIntensityMap()
	{	return &intensityMap_;};
	std::map<double, double>* getIntensityHist()
	{	return &intensityHist_;};
	double getMean()
	{	return mean_;};
};

} // ns

#endif

