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
// $Maintainer: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERHIRES_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <map>

#define DEBUG_PEAK_PICKING
#undef DEBUG_PEAK_PICKING
//#undef DEBUG_DECONV
namespace OpenMS
{
  /**
		 @brief This class makes nothing.
		 
		 @htmlinclude OpenMS_PeakPickerCWT.parameters
		 
		 @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerHiRes
		: public DefaultParamHandler,
			public ProgressLogger
  {
	public:
    /// Constructor
    PeakPickerHiRes();
		
    /// Destructor
    virtual ~PeakPickerHiRes();
		
    /** 
				@brief pick
    */
    template <typename PeakType>
    void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)
    {
			// double SNT_THRESHOLD = 2.0;
			
			// signal2noise estimation
			SignalToNoiseEstimatorMedian<MSSpectrum<PeakType> > snt;
			snt.init(input);
			
			// find local maxima in raw data
			for (unsigned int i = 1; i < input.size()-1; ++i)
				{
					double centralPeakMZ = input[i].getMZ(), centralPeakInt = input[i].getIntensity();
					double lNeighbourMZ = input[i-1].getMZ(), lNeighbourInt = input[i-1].getIntensity();
					double rNeighbourMZ = input[i+1].getMZ(), rNeighbourInt = input[i+1].getIntensity();
					
					// MZ spacing sanity checks
					double leftToCentral = std::fabs(centralPeakMZ-lNeighbourMZ);
					double centralToRight = std::fabs(rNeighbourMZ-centralPeakMZ);
					double minSpacing = (leftToCentral < centralToRight)? leftToCentral : centralToRight;
				
					// look for peak cores meeting MZ and intensity/SNT criteria
					if (snt.getSignalToNoise(input[i]) >= signalToNoise_
							&& leftToCentral < 1.5*minSpacing
							&& centralPeakInt > lNeighbourInt
							&& snt.getSignalToNoise(input[i-1]) >= signalToNoise_
							&& centralToRight < 1.5*minSpacing
							&& centralPeakInt > rNeighbourInt
							&& snt.getSignalToNoise(input[i+1]) >= signalToNoise_)
						{
							// special case: if a peak core is surrounded by more intense
							// satellite peaks (indicates oscillation rather than
							// real peaks) -> remove							
							if ((i > 1
									 && std::fabs(lNeighbourMZ - input[i-2].getMZ()) < 1.5*minSpacing
									 && lNeighbourInt < input[i-2].getIntensity()
									 && snt.getSignalToNoise(input[i-2]) >= signalToNoise_)
									&&
									((i+2) < input.size()
									 && std::fabs(input[i+2].getMZ() - rNeighbourMZ) < 1.5*minSpacing
									 && rNeighbourInt < input[i+2].getIntensity()
									 && snt.getSignalToNoise(input[i+2]) >= signalToNoise_)
									)
								{
									++i;
									continue;
								}
							
							
							std::map<double, double> peakRawData;
							
							peakRawData[centralPeakMZ] = centralPeakInt;
							peakRawData[lNeighbourMZ] = lNeighbourInt;
							peakRawData[rNeighbourMZ] = rNeighbourInt;
							
							
							// peak core found, now extend it
							// to the left
							unsigned int k = 2;
							
							while ((i-k+1) > 0
										 && std::fabs(input[i-k].getMZ() - peakRawData.begin()->first) < 1.5*minSpacing
										 && input[i-k].getIntensity() < peakRawData.begin()->second
										 && snt.getSignalToNoise(input[i-k]) >= signalToNoise_)
								{
									peakRawData[input[i-k].getMZ()] = input[i-k].getIntensity();
									++k;
								}
							
							// to the right
							k = 2;
							while ((i+k) < input.size()
										 && std::fabs(input[i+k].getMZ() - peakRawData.rbegin()->first) < 1.5*minSpacing
										 && input[i+k].getIntensity() < peakRawData.rbegin()->second
										 && snt.getSignalToNoise(input[i+k]) >= signalToNoise_)
								{
									peakRawData[input[i+k].getMZ()] = input[i+k].getIntensity();
									++k;
								}
						  
							
							// output all raw data points selected for one peak
							// TODO: #ifdef DEBUG_ ...
							/**
								 for (std::map<double, double>::const_iterator mIt = peakRawData.begin(); mIt != peakRawData.end(); ++mIt) {
								 PeakType peak;
								 
								 peak.setMZ(mIt->first);
								 peak.setIntensity(mIt->second);
								 output.push_back(peak);
								 
								 std::cout << mIt->first << " " << mIt->second << " snt: " << std::endl;
								 
								 }
								 std::cout << "--------------------" << std::endl;
							*/
							
							
							const size_t numRawPoints = peakRawData.size();
							
							std::vector<double> rawMZvalues;
							std::vector<double> rawINTvalues;
							
							for (std::map<double, double>::const_iterator mIt = peakRawData.begin(); mIt != peakRawData.end(); ++mIt)
								{
									rawMZvalues.push_back(mIt->first);
									rawINTvalues.push_back(mIt->second);
								}
							
							// setup gsl splines
							gsl_interp_accel *splineAcc = gsl_interp_accel_alloc();
							gsl_interp_accel *deriv1Acc = gsl_interp_accel_alloc();
							gsl_spline *peakSpline = gsl_spline_alloc(gsl_interp_cspline, numRawPoints);
							gsl_spline_init(peakSpline, &(*rawMZvalues.begin()), &(*rawINTvalues.begin()), numRawPoints);
							
							
							// calculate maximum by evaluating the spline's 1st derivative
							// (bisection method)
							double maxMZ = centralPeakMZ, maxInt = centralPeakInt;
							double threshold = 0.000001;
							double l = lNeighbourMZ;
							double r = rNeighbourMZ;
							
							bool lSign = 1;
							double eps = std::numeric_limits<double>::epsilon();

							
							// bisection
							do
								{
									double mid = (l + r) / 2;
									
									double fMid = gsl_spline_eval_deriv(peakSpline, mid, deriv1Acc);
									if (!(std::fabs(fMid) > eps))
										break;

									bool mSign = (fMid < 0.0)? 0 : 1;
									
									if (lSign ^ mSign)
										{
											r = mid;
										}
									else
										{
											l = mid;
										}
									
									// TODO: #ifdef DEBUG_ ...
									// PeakType peak;
									// peak.setMZ(mid);
									// peak.setIntensity(gsl_spline_eval(peakSpline, mid, splineAcc));
									// output.push_back(peak);
									
								} while ( std::fabs(l - r) > threshold );
							
							// sanity check?
							maxMZ = (l + r) / 2;
							maxInt = gsl_spline_eval(peakSpline, maxMZ, splineAcc);
							
							// save picked pick into output spectrum
							PeakType peak;
							peak.setMZ(maxMZ);
							peak.setIntensity(maxInt);
							output.push_back(peak);
							
							// free allocated gsl memory
							gsl_spline_free(peakSpline);
							gsl_interp_accel_free(splineAcc);
							gsl_interp_accel_free(deriv1Acc);
							
							// jump over raw data points that have been considered already
							i = i + k - 1;
						}
				}
			
			return ;
    }
		
		
    /** 
				@brief 
    */
    template <typename PeakType>
    void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)
    {
			for (unsigned int scanIdx = 0; scanIdx != input.size(); ++scanIdx)
				{
					if (input[scanIdx].getMSLevel() == 1)
						{
							MSSpectrum<PeakType> tmpSpectrum;
							
							pick(input[scanIdx], tmpSpectrum);
							
							tmpSpectrum.setRT(input[scanIdx].getRT());
							output.push_back(tmpSpectrum);
						}
					// post a warning if MS2 scans are detected?
				}
			
			return ;
    }
		
	protected:
		// signal-to-noise parameter
		double signalToNoise_;
		
		
		//docu in base class
		void updateMembers_();
		
	}; // end PeakPickerHiRes
	
}// namespace OpenMS

#endif
