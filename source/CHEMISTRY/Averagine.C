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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/Averagine.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>

	
namespace OpenMS
{

double Averagine::gay_constants[6][5] 
	= {	{M01,M02,M03,M04,M05},
			{M11,M12,M13,M14,M15}, 
			{M21,M22,M23,M24,M25}, 
			{M31,M32,M33,M34,M35},
			{M41,M42,M43,M44,M45},
			{M51,M52,M53,M54,M55} };


std::vector<RawDataPoint1D> Averagine::getModel (const double mz_pos, const unsigned int charge, 
	std::pair<double, double>* pattern_extend, const double cut_off) throw ()
{
	const double mass = mz_pos*charge - charge*EXACT_NEUTRON_MASS;
	std::vector<double> gay_mass (5);
	gay_mass[4] = 1.;
	gay_mass[3] = (mass-800.)/2200.;
	gay_mass[2] = gay_mass[3]*gay_mass[3];
	gay_mass[1] = gay_mass[2]*gay_mass[3];
	gay_mass[0] = gay_mass[1]*gay_mass[1];
	
	std::vector<RawDataPoint1D> res (6);	

	double tmp;
	for (unsigned int p=0; p<res.size(); ++p) //maybe we should think about delooping and hard coding for a substantial speed up
	{
		tmp=0;
		for (unsigned int i=0; i<5; ++i)
		{ 
			tmp += gay_constants[p][i]*gay_mass[i]; 
		};
		
		res[p].setIntensity (tmp);
		res[p].setMZ (mz_pos+(NEUTRON_MASS*p/(double)charge));	
	};

	if (pattern_extend != NULL)
	{
		std::vector<RawDataPoint1D>::iterator below = res.begin()+3;
		while (below->getIntensity() > cut_off && below != res.end())
		{
			++below;
		};
		
		pattern_extend->first = mz_pos - 0.25*NEUTRON_MASS/(double)charge;
	
		if (below == res.end()) //i.e. all peaks are significant
		{
			pattern_extend->second = mz_pos + 5.5*NEUTRON_MASS/(double)charge;
		}
		else
		{
			pattern_extend->second = (below-1)->getMZ() + 0.5*NEUTRON_MASS/(double)charge;
		};
	};

	return(res);
}

} // namespace OpenMS

