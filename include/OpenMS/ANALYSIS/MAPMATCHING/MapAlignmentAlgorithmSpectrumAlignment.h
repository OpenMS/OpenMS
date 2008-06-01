// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Vipul Patel $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on spectrum similarity (dynamic programming). 		
	*/
	class MapAlignmentAlgorithmSpectrumAlignment
	 : public MapAlignmentAlgorithm, 
	 	public ProgressLogger
	{
		public:
			/// Default constructor
			MapAlignmentAlgorithmSpectrumAlignment();

			/// Destructor
			virtual ~MapAlignmentAlgorithmSpectrumAlignment();

			//Docu in base class
		virtual void alignPeakMaps(std::vector< MSExperiment<> >&);
			
			
			///Creates a new instance of this class (for Factory)
			static MapAlignmentAlgorithm* create()
			{
				return new MapAlignmentAlgorithmSpectrumAlignment();
			}
			
			///Returns the product name (for the Factory)
			static String getProductName()
			{
				return "spectrum_alignment";
			}
			
		private:
			///Copy constructor is not implemented -> private
			MapAlignmentAlgorithmSpectrumAlignment(const MapAlignmentAlgorithmSpectrumAlignment& );
			///Assignment operator is not implemented -> private
			MapAlignmentAlgorithmSpectrumAlignment& operator=(const MapAlignmentAlgorithmSpectrumAlignment& );
			void prepareAlign_(const std::vector< MSSpectrum<>* >&, MSExperiment<>& );
			void msFilter_(MSExperiment<>& peakmap,std::vector<MSSpectrum<>* >& spectrum_pointer_container);
			void fourierActivation_(std::vector<MSSpectrum<>* >& spectrum_pointer_container);
			void transform_(MSSpectrum<> & spec);
			bool insideBand_(UInt i,Int j,UInt n,UInt m,Int k_);
			void ordering_(std::vector<int>& x,std::vector<double>& y );
			void calculateSpline_(std::vector<int>& x,std::vector<double>& y, std::vector<MSSpectrum<>* >& aligned,UInt begin, UInt end);
			Int bestk_(const std::vector<MSSpectrum<>* >& pattern, std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation, UInt xbegin,UInt xend, UInt ybegin, UInt yend);
			Real scoreCalculation_(UInt i,Int j, UInt patternbegin, UInt alignbegin ,const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation);
			Real scoring_(const MSSpectrum<>& a, MSSpectrum<>& b);
			void affineGapalign_(UInt xbegin, UInt ybegin, UInt xend,UInt yend, const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::vector<int> xcordinate, std::vector<double>ycordinate);
			///gap-cost
			Int gap_;
			///extension cost
			Int e_;
			///
			PeakSpectrumCompareFunctor* c1_;
			///threshold
			Real cutoffScore_;


		

	void updateMembers_();


	};
	 
	

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H



