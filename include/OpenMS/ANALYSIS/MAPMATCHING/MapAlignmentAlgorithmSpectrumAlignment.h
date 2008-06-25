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
		
	  @ref MapAlignmentAlgorithmSpectrumAlignment_Parameters are explained on a separate page.  
	  
	  @todo fix warnings (Vipul)
	  @todo fill transformations with pairs as soon as this is implemented in TransformationDescription (Vipul)
	  
		@experimental This algorithm is work in progress and might change.

		@ingroup MapAlignment
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
			virtual void alignPeakMaps(std::vector< MSExperiment<> >&, std::vector<TransformationDescription>&);
			
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
			
			/**
				@brief A function to prepare the sequence for the alignment. It calls intern the main function for the alignment.
		   
				This function takes two arguments. These argument types are two MSExperiments. 
				The first argument should have been filtered, so that only the type of MSLevel 1 exists in the Sequence. 
				The second argument doesn't have to fulfill this restriction. It's going to be filtered automatically. 
				With these two arguments a precalculation is done to find some corresponding data points(maximum 4) for building alignment blocks.
				After the alignment a retransformation is done, the new Retention Times appear in the original data.      
		  
			  The parameters are MSExperiments.
			  @param pattern is the template MSExperiment.
			  @param aligned is the sequence which has to be aligned(also MSExperiment).
			  @see MapAlignmentAlgorithmSpecturmAlignment()
			*/	
			
			void prepareAlign_(const std::vector< MSSpectrum<>* >&, MSExperiment<>& );
			/**
			  @brief filtered the MSLevel to gain only MSLevel 1 
			
			 	The alignment works only on MSLevel 1 data, so a filter has to be run.
			 
				@param peakmap MSExperiment is the sequence which has to be filtered
		 		@param spectrum_pointer_container std::vector<MSSpectrum<>* > output container, where pointers of the MSSpectrum are saved(only with MSLevel 1 type)
			*/		
		
			void msFilter_(MSExperiment<>& peakmap,std::vector<MSSpectrum<>* >& spectrum_pointer_container);
			/**
				@brief does the transformation if the Discrete Cosines Fourier Transformation is selected.
		
			 	Call internally the function transform, only if the comparison score function Fourier is selected.
				@param spectrum_pointer_container is the sequence which has to be transform
			*/			
			
			void fourierActivation_(std::vector<MSSpectrum<>* >& spectrum_pointer_container);
			/**
		  	@brief calculate the Discrete Cosines Fourier Transformation.
		     				
		 		This Function transforms a given MSSpectrum to a Discrete Cosines Fourier Transformation. It stores only the part of the cosines of the FFT in
		 		the MetaDataArray which is a container from the MSSpectrum. Only call this function, if you are sure there is no other 				transformation done earlier over the same MSSpectrum, because it isn't checked if there already exists a transformation.
		        		
		   	@param spec  MSSpectrum 
		  */			
			
			void transform_(MSSpectrum<> & spec);
			/**
				@brief function for the test if cell i,j of the grid is inside the band

				The function returns true if the cell underlie these conditions:
				-k<=i-j<=k+n-m
				else retun false.
				@param i Int coordinate i
				@param j Int coordinate j
				@param n UInt  size of column
				@param m UInt  size of row
				@param k_ Int  size of k_
			*/			
			bool insideBand_(UInt i,Int j,UInt n,UInt m,Int k_);
			
			//void ordering_(std::vector<int>& x,std::vector<double>& y );
			/**
		    @brief function to calculate a cubicspline to interpolate the Retention time

			 	calculateSpline_(cubic spline) is needed to interpolate the Retention time for the whole length of the sequence.	
			 	The data points are the points in which a match appeared. To get the rest of the Retention times a spline is necessary.
			 	
			 	@param x std::vector<int> which contain x cordinates
			 	@param y  std::vector<double> which contain the retentiontimes
			 	@param aligned MSExperiment the aligned sequence
			 	@param begin UInt begin of the alignment in the aligned sequence 
			 	@param end UInt end of the alignment in the aligned sequence 
			 	@see MapAlignmentAlgorithmSpectrumAlignment()
			*/			
			void calculateSpline_(std::vector<int>& x,std::vector<double>& y, std::vector<MSSpectrum<>* >& aligned,UInt begin, UInt end);
		  /**
			 	@brief calculate the size of the band for the alignment for two given Sequence
			   		
     		This function calculates the size of the band for the alignment. It takes three samples from the aligned sequence and tries to 
     		find the highscore pairs(matching against the template sequence). The highscore pair with the worst distance is to be chosen as the size of k. 
     		
     		@param pattern const std::vector<T* > vector of pointers of the template sequence
     		@param aligned std::vector<T* > vector of pointers of the aligned sequence
     		@param buffer std::map<UInt, std::map<UInt,Real> >  holds the calculated score of index i,j.
     		@param flag1 bool flag indicate the order of the matrix   		
     		@param xbegin UInt indicate the beginning of the template sequence
     		@param xend UInt indicate the end of the template sequence
     		@param ybegin UInt indicate the beginning of the aligned sequence
     		@param yend UInt indicate the end of the aligned sequence
     	   	
    	*/			
			Int bestk_(const std::vector<MSSpectrum<>* >& pattern, std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation, UInt xbegin,UInt xend, UInt ybegin, UInt yend);
			/**
				@brief calculate the score of two given MSSpectrums calls intern scoring_

				This function calculates the score from two MSSpectrums. These two MSSpectrums are chosen by the coordinates i,j. 
				I,j indicate the index in the matrix. To find the right index on the sequence, each beginning is also given to the function.
				A flag indicates the labeling of the axes. The buffermatrix stores the result of the scoring. If the band expands only a lookup of known scores is done.

				@param i UInt is a index from the matrix.
				@param j Int is a index from the matrix.
				@param patternbegin UInt indicate the beginning of the template sequence
				@param alignbegin UInt indicate the beginning of the aligned sequence
				@param pattern const std::vector<T* > vector of pointers of the template sequence
				@param aligned std::vector<T* > vector of pointers of the aligned sequence
				@param buffer std::map<UInt, std::map<UInt,Real> >  holds the calculated score of index i,j.
				@param flag1 bool flag indicate the order of the matrix
				
			*/
			Real scoreCalculation_(UInt i,Int j, UInt patternbegin, UInt alignbegin ,const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::map<UInt, std::map<UInt,Real> > & buffer,bool column_row_orientation);
			/**
				@brief return the score of two given MSSpectrums by calling the scorefunction
		 		
		 		@param a is a MSSpectrum.
			 	@param b is a MSSpectrum
			*/			
			Real scoring_(const MSSpectrum<>& a, MSSpectrum<>& b);
			/**
				@brief affine gap cost Alignment

				This Alignment is based on the Needleman Wunsch Algorithm. 
				To improve the time complexity a banded version was implemented, known as k - alignment. 
				To save some space, the alignment is going to be calculated by position xbegin to xend of one sequence and ybegin 
				and yend by another given sequence. The result of the alignment is stored in the second argument. 
				The first sequence is used as a template for the alignment.

				 	 
			 	@param xbegin UInt cordinate for the beginning of the template sequence.
			 	@param ybegin UInt cordinate for the beginning of the aligend sequence .
			 	@param xend UInt cordinate for the end of the template sequence.
			 	@param yend UInt cordinate for the end of the aligend sequence.
			 	@param pattern is the template MSExperiment.
			 	@param aligned is to be align MSExperiment.
			 	@see MapAlignmentAlgorithmSpecturmAlignment()
		 */			
			void affineGapalign_(UInt xbegin, UInt ybegin, UInt xend,UInt yend, const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::vector<int>& xcordinate, std::vector<double>&ycordinate);
			
			///Represent the gap cost for opening o closing a gap in the alignment
			Int gap_;
			///This is the extension cost after a gap ist open
			Int e_;
			///Pointer holds the scoringfunction, which can be selected
			PeakSpectrumCompareFunctor* c1_;
			///This is the minimal score to be count as a mismatch(range 0.0 - 1.0)
			Real cutoffScore_;
			

			void updateMembers_();

	};
	 
	

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H



