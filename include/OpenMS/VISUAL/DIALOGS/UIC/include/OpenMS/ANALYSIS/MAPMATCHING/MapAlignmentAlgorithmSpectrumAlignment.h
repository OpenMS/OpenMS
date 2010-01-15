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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fft_real.h>

namespace OpenMS
{
	/**
		@brief A map alignment algorithm based on spectrum similarity (dynamic programming). 		
		
	  @htmlinclude OpenMS_MapAlignmentAlgorithmSpectrumAlignment.parameters
	  
		@experimental This algorithm is work in progress and might change.
		
	
		@ingroup MapAlignment
	*/
	class OPENMS_DLLAPI MapAlignmentAlgorithmSpectrumAlignment
	 : public MapAlignmentAlgorithm
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
				@brief innerclass necessary for using the sort algo. 		
							
				Defines serveral overloaded operator() for using the std::sort 	algorithm 
				
				@experimental This algorithm is work in progress and might change.
		
			*/
			class OPENMS_DLLAPI Compare
			{
				protected:
					bool flag;
				
				public:
					
					/// Default constructor with an order flag
					Compare(bool b=false)
					 : flag(b)
					{
					}
					///overloaded operator() for comparing maps of maps std::pair<std::pair<Int,Real>,Real>. If flag is false the second argument of the outer map is selected. The output is an ascending order. If the order flag is true, the first argument of the inner class is selected to get a descending oder.
					inline bool operator()(const std::pair<std::pair<Int,Real>,Real>& c1, const std::pair<std::pair<Int,Real>,Real>& c2 )
					{
						if(!flag)
						{
							return c1.second > c2.second;
						}
						else
						{
							return (c1.first).first < (c2.first).first;
						}
					}
					///overloaded operator() for comparing pairs of Real, Real std::pair<Real,Real>. If the order flag is false, an ascending order are returned else a descending. The comparison is done by the first argument of the map.
					inline bool operator()(const std::pair<Real,Real> & c1, const std::pair<Real,Real>& c2)
					{
						if(!flag)
						{
							return c1.first > c2.first;
						}
						else
						{
							return c1.first < c2.first;
						}
					}
			};
			
			/**
				@brief A function to prepare the sequence for the alignment. It calls intern the main function for the alignment.
		   
				This function takes two arguments. These argument types are two MSExperiments. 
				The first argument should have been filtered, so that only the type of MSLevel 1 exists in the Sequence. 
				The second argument doesn't have to fulfill this restriction. It's going to be filtered automatically. 
				With these two arguments a precalculation is done to find some corresponding data points(maximum 4) for building alignment blocks.
				After the alignment a retransformation is done, the new Retention Times appear in the original data.      
		  
			  The parameters are MSExperiments.
				@param pattern template map.
				@param aligned map which has to be aligned.
				@param transformation container for rebuilding the alignment only by specific data-points
			*/	
			void prepareAlign_(const std::vector< MSSpectrum<>* >& pattern, MSExperiment<>& aligned,std::vector<TransformationDescription>& transformation );
			
			/**
			  @brief filtered the MSLevel to gain only MSLevel 1 
			
			 	The alignment works only on MSLevel 1 data, so a filter has to be run.
			 
				@param peakmap map which has to be filtered
		 		@param spectrum_pointer_container output container, where pointers of the MSSpectrum are saved(only with MS level 1)
			
				@exception Exception::IllegalArgument is thrown if no spectra are contained in @p peakmap
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
		 		the FloatDataArray which is a container from the MSSpectrum. Only call this function, if you are sure there is no other 				transformation done earlier over the same MSSpectrum, because it isn't checked if there already exists a transformation.
		  */			
			void transform_(MSSpectrum<> & spec);

			/**
				@brief function for the test if cell i,j of the grid is inside the band

				The function returns true if the cell underlie these conditions:
				-k<=i-j<=k+n-m
				else retun false.
				@param i coordinate i
				@param j coordinate j
				@param n size of column
				@param m size of row
				@param k_ size of k_
			*/			
			bool insideBand_(Size i,Size j,Size n,Size m,Int k_);
			
			/**
		    @brief function to calculate a cubicspline to interpolate the Retention time

			 	calculateSpline_(cubic spline) is needed to interpolate the Retention time for the whole length of the sequence.	
			 	The data points are the points in which a match appeared. To get the rest of the Retention times a spline is necessary. The result of the spline is saved int the TransformationDescription Container.
			 	
			 	@param x the coordinate
			 	@param y the retention times
			 	@param aligned map
			 	@param begin of the alignment in the aligned sequence 
			 	@param end of the alignment in the aligned sequence 
				@param transformation specific data points to recalulcate the spline
			*/			
			void calculateSpline_(std::vector<int>& x,std::vector<Real>& y, std::vector<MSSpectrum<>* >& aligned,Size begin, Size end,std::vector<TransformationDescription>& transformation);
		  
		  /**
			 	@brief calculate the size of the band for the alignment for two given Sequence
			   		
     		This function calculates the size of the band for the alignment. It takes three samples from the aligned sequence and tries to 
     		find the highscore pairs(matching against the template sequence). The highscore pair with the worst distance is to be chosen as the size of k. 
     		
				@param pattern vector of pointers of the template sequence
				@param aligned vector of pointers of the aligned sequence
     		@param buffer holds the calculated score of index i,j.
     		@param column_row_orientation indicate the order of the matrix 	
     		@param xbegin indicate the beginning of the template sequence
     		@param xend indicate the end of the template sequence
     		@param ybegin indicate the beginning of the aligned sequence
     		@param yend indicate the end of the aligned sequence
    	*/			
			Int bestk_(const std::vector<MSSpectrum<>* >& pattern, std::vector<MSSpectrum<>* >& aligned,std::map<Size, std::map<Size,Real> > & buffer,bool column_row_orientation, Size xbegin,Size xend, Size ybegin, Size yend);
			
			/**
				@brief calculate the score of two given MSSpectrums calls intern scoring_

				This function calculates the score from two MSSpectrums. These two MSSpectrums are chosen by the coordinates i,j. 
				I,j indicate the index in the matrix. To find the right index on the sequence, each beginning is also given to the function.
				A flag indicates the labeling of the axes. The buffermatrix stores the result of the scoring. If the band expands only a lookup of known scores is done.

				@param i is a index from the matrix.
				@param j is a index from the matrix.
				@param patternbegin indicate the beginning of the template sequence
				@param alignbegin  indicate the beginning of the aligned sequence
				@param pattern vector of pointers of the template sequence
				@param aligned vector of pointers of the aligned sequence
				@param buffer  holds the calculated score of index i,j.
				@param column_row_orientation indicate the order of the matrix
			*/
			Real scoreCalculation_(Size i,Size j, Size patternbegin, Size alignbegin ,const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::map<Size, std::map<Size,Real> > & buffer,bool column_row_orientation);
			
			/**
				@brief return the score of two given MSSpectrums by calling the scorefunction
			*/			
			Real scoring_(const MSSpectrum<>& a, MSSpectrum<>& b);
			
			/**
				@brief affine gap cost Alignment

				This Alignment is based on the Needleman Wunsch Algorithm. 
				To improve the time complexity a banded version was implemented, known as k - alignment. 
				To save some space, the alignment is going to be calculated by position xbegin to xend of one sequence and ybegin 
				and yend by another given sequence. The result of the alignment is stored in the second argument. 
				The first sequence is used as a template for the alignment.

				 	 
			 	@param xbegin cordinate for the beginning of the template sequence.
			 	@param ybegin cordinate for the beginning of the aligend sequence .
			 	@param xend cordinate for the end of the template sequence.
			 	@param yend cordinate for the end of the aligend sequence.
				@param pattern template map.
				@param aligned map to be aligned.
			 	@param xcoordinate save the postion of ankerpoints
   			@param ycoordinate save the retentiontimes of an ankerpoints
   			@param xcoordinatepattern save the reference position of the ankerpoints from the pattern
			
				@exception Exception::OutOfRange if a out of bound appear @p pattern or @p aligned
			*/
			void affineGapalign_(Size xbegin, Size ybegin, Size xend,Size yend, const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned,std::vector<int>& xcoordinate, std::vector<Real>&ycoordinate, std::vector<int>& xcoordinatepattern);
			
			/**
				@brief  preparation function of data points to construct later the  spline function.

				This function reduced the amount of data values for the next step. The reduction is done by using a number of buckets, where the data points a selected. 
				Within the buckets, only defined number a selected, to be written back as a data point. 
				The selection within the buckets is done by scoring.
					
				@param pattern template map.
				@param aligned map to be aligned.
				@param xcoordinate save the position of anchor points
				@param ycoordinate  save the retention times of an anchor points
				@param xcoordinatepattern save the reference position of the anchor points from the pattern
			*/
			void bucketFilter_(const std::vector<MSSpectrum<>* >& pattern,std::vector<MSSpectrum<>* >& aligned,std::vector<Int> & xcoordinate, std::vector<Real> & ycoordinate, std::vector<Int>&xcoordinatepattern);
			
			/**
				@brief Creates files for the debugging

				This function is only active if the debugflag ist true. The debugfileCreator creates following files debugtraceback.txt(gnuplotScript), debugscoreheatmap.r and debugRscript. Debugscoreheatmap.r contains the scores of the Spectra to each other from the alignment and also the traceback. DebugRscript is the R script which reads those data. So both files are only working under R. Start R and type main(location of debugscoreheatmap.r). The output will be a heatmap of each sub-alignment. Debugtraceback.txt shows the way of the Traceback by using gnuplot.

				@param pattern template map.
				@param aligned map to be aligned.
			*/
			void debugFileCreator_(const std::vector<MSSpectrum<>* >& pattern,  std::vector<MSSpectrum<>* >& aligned);
			
			/**
				@brief Delete entries of the FloatDataArray which was made from CompareFouriertransform 

				This function erase the entries with was done by the CompareFouriertransform function.
			*/
			void eraseFloatDataArrayEntry_(std::vector<MSSpectrum<>* >& spectrum_pointer_container);
			
			/**
				@brief Rounding the score of two spectras, only necessary for debugging
				
				This function rounded the score of two spectras.This is necessary for some function in the Debug-Mode
			*/
			void debugscoreDistributionCalculation_(Real score);
			///Represent the gap cost for opening o closing a gap in the alignment
			Real gap_;
			///Extension cost after a gap ist open
			Real e_;
			///Pointer holds the scoringfunction, which can be selected
			PeakSpectrumCompareFunctor* c1_;
			///This is the minimal score to be count as a mismatch(range 0.0 - 1.0)
			Real cutoffScore_;
			///Defines the size of one bucket 
			Size bucketsize_;
			///Defines the amount of ankerpoints which are selected within one bucket. 
			Size anchorPoints_;
			///Debug mode flag default: False
			bool debug_;
			///Represent the cost of a mismath in the alignment
			Real mismatchscore_;
			///This is the minimum score for counting as a match(1-cutoffScore_)
			Real threshold_;
			///Container holding the score of the matchmatrix and also the insertmatrix
			std::vector<std::vector<Real> >debugmatrix_;
			///Container holding the only the score of Spectrums
			std::vector<std::vector<Real> >debugscorematrix_;
			///Container holding the path of the traceback
			std::vector<std::pair<Real,Real> >debugtraceback_;
			///Container holding the score of each cell(matchmatrix,insertmatrix, traceback)
			std::vector<Real> scoredistribution_;//save the cell i, j , matchscore, insertscore, traceback
			//docu in base class
			void updateMembers_();
	};
	 
	

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTALGORITHMSPECTRUMALIGNMENT_H




