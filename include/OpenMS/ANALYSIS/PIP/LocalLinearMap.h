// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Alexandra Scherbart $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_PIP_LOCALLINEARMAP_H
#define OPENMS_ANALYSIS_PIP_LOCALLINEARMAP_H
	
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
	
namespace OpenMS
{

	/**
		@brief Trained Local Linear %Map (LLM) model for peak intensity prediction
		
		This class offers a model for predictions of peptide peak heights 
		(referred to as intensities) by a Local Linear %Map (LLM) model and
		is the basis of PeakIntensityPredictor.
		
		A general introduction to the Peak Intensity Predictor (PIP)
		can be found in the <A HREF="tutorial_pip.html">PIP Tutorial</A>.
		
		The model trained needs two files for storing the position of the
		codebook vectors and the linear mappings (codebooks.data, linearMapping.data)
		This is the default model used by PeakIntensityPredictor.
  */	
	class LocalLinearMap
	{
	
		public:	  
		
			/**
	   		@brief Define parameters needed by the Local Linear %Map (LLM) model
	    
	  		Parameters xdim and ydim define the size of the two dimensional
	  		grid structure. Parameter radius gives the width of the Gaussian
	  		neighborhood function.
  		*/	
			struct LLMParam 
			{
				UInt xdim; /**< size of first coordinate */
				UInt ydim; /**< size of second coordinate */
				DoubleReal radius; /**< width of Gaussian neighborhood function */
			};
		
			/** @name Constructors and Destructors
			*/
			//@{
			/// default constructor
			LocalLinearMap();
			/**
			 @brief constructor with default model load
			
				Load constructor with default model
				@param codefile  read file with positions of prototypes (codebooks)
				@param Afile read file corresponding to linear mappings
				@param xdim size of first coordinate of grid structure
				@param ydim size of second coordinate of grid structure
				@param radius width of Gaussian neighborhood
				@see LLMParam
			 */
			LocalLinearMap(String codefile, String Afile, UInt xdim, UInt ydim, DoubleReal radius);
			/// destructor
			virtual ~LocalLinearMap();
			//@}
	
			///return parameters of the LocalLinearMap model
			const LLMParam& getLLMParam() const;
			///return position of the codebook vectors (18-dim)
			const Matrix<DoubleReal>& getCodebooks() const;
			///return linear mappings of the codebooks
			const Matrix<DoubleReal>& getMatrixA() const;
			///return linear bias
			const std::vector< DoubleReal >& getVectorWout() const;
			///return coordinates of codebook vectors on the 2-d grid
			const Matrix<UInt>& getCord() const;
			///calculate and return normalized amino acid index variables from string representation of peptide
			void normalizeVector(std::vector<DoubleReal>& aaIndexVariables);
			///calculate neighborhood function based on distance of prototypes to winner prototype on two-dimensional grid structure and neighborhood width.
			std::vector<DoubleReal> neigh(const Matrix<UInt>& cord, UInt win, DoubleReal radius);
			
		private:
			
			LLMParam param_;									///<parameters of the model
			Matrix<DoubleReal> code_;					///<codebook vectors
			Matrix<DoubleReal> A_;						///<linear mappings
			std::vector< DoubleReal > wout_;	///<linear bias
			Matrix<UInt> cord_;								///<coordinates of codebooks on grid

			/// needed to store prototype coordinates
			Matrix<UInt> genCord_(UInt xdim, UInt ydim);
			///calculate distance between two prototypes
			DoubleReal dist_(const Matrix<UInt>& u, const Matrix<UInt>& v, UInt a, UInt b);
			
			///Copy constructor not implemented => private
			LocalLinearMap(LocalLinearMap& rhs);
			///Assignment operator not implemented => private
			LocalLinearMap& operator = (const LocalLinearMap& llm);
		
		};
	

}
#endif 
