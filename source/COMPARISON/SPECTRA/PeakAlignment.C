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
// $Maintainer: Mathias Walzer $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>

using namespace std;

namespace OpenMS
{
	PeakAlignment::PeakAlignment()
	  : PeakSpectrumCompareFunctor()
	{
			defaults_.setValue("epsilon", 0.2, "defines the absolut error of the mass spectrometer");
			defaults_.setValue("normalized", 1, "is set 1 if the similarity-measurement is normalized to the range [0,1]");
			defaults_.setValue("precursor_mass_tolerance", 3.0, "Mass tolerance of the precursor peak, defines the distance of two PrecursorPeaks for which they are supposed to be from different peptides");
			defaultsToParam_();
	}
	
	PeakAlignment::PeakAlignment(const PeakAlignment& source)
	  : PeakSpectrumCompareFunctor(source)
	{
	}
	
	PeakAlignment::~PeakAlignment()
	{
	}
	
	PeakAlignment& PeakAlignment::operator = (const PeakAlignment& source)
	{
		if (this != &source)
		{
	  		PeakSpectrumCompareFunctor::operator = (source);
		}
	  	return *this;
	}

	double PeakAlignment::operator () (const PeakSpectrum& spec) const
	{
		return operator () (spec, spec);
	}
	
	double PeakAlignment::operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const
	{
		
		// shortcut similarity calculation by comparing PrecursorPeaks (PrecursorPeaks more than delta away from each other are supposed to be from another peptide)
	  	const double delta = (double)param_.getValue("precursor_mass_tolerance");
	  	if(fabs(spec1.getPrecursorPeak().getPosition()[0]-spec2.getPrecursorPeak().getPosition()[0])>delta)
	  	{
	  		return 0;
	  	}
	  	
		const double epsilon = (double)param_.getValue("epsilon");
		
		//TODO gapcost dependence on distance ?
		const double gap = (double)param_.getValue("epsilon");

		//initialize alignment matrix with 0 in (0,0) and a multiple of gapcost in the first row/col matrix(row,col,values)
		Matrix<double> matrix(spec1.getContainer().size()+1, spec2.getContainer().size()+1,0);
		for (UInt i = 1; i < matrix.rows(); i++)
		{
			matrix.setValue(i,0,-gap*i);
		}
		for (UInt i = 1; i < matrix.cols(); i++)
		{
			matrix.setValue(0,i,-gap*i);
		}		
		
		//get sigma - the standard deviation (sqrt of variance)
		double mid(0);
		for (UInt i = 0; i < spec1.getContainer().size(); ++i)
		{
			for (UInt j = 0; j < spec2.getContainer().size(); ++j)
			{
				double pos1(spec1.getContainer()[i].getMZ()), pos2(spec2.getContainer()[j].getMZ());
				mid += fabs(pos1 - pos2);
			}
		}
		// average peak distance
		mid /= (spec1.getContainer().size()*spec2.getContainer().size());
			
			#ifdef OPENMS_DEBUG 	
			cout << "average peak distance " << mid << endl; 
			#endif
				
				
		double var(0);
		for (UInt i = 0; i < spec1.getContainer().size(); ++i)
		{
			for (UInt j = 0; j < spec2.getContainer().size(); ++j)
			{
				double pos1(spec1.getContainer()[i].getMZ()), pos2(spec2.getContainer()[j].getMZ());
				var += (fabs(pos1 - pos2) - mid)*(fabs(pos1 - pos2) - mid);
			}
		}	
		// peak distance variance
		var = var/(spec1.getContainer().size()*spec2.getContainer().size());
				
			#ifdef OPENMS_DEBUG 
			cout << "peak distance variance " << var << endl;
			#endif

		//only in case of only two equal peaks in the spectra sigma is 0 	


		const double sigma((var==0)?FLT_MIN:sqrt(var));
				
			#ifdef OPENMS_DEBUG 
			cout << "peak standard deviation " << sigma << endl;
			#endif
				
		//fill alignment matrix
		for (UInt i = 1; i < spec1.getContainer().size()+1; ++i)
		{
			for (UInt j = 1; j < spec2.getContainer().size()+1; ++j)
			{
				double pos1(spec1.getContainer()[i-1].getMZ()), pos2(spec2.getContainer()[j-1].getMZ());
				//only if peaks are in reasonable proximity alignment is considered else only gaps
				if (fabs(pos1 - pos2) <= epsilon)
				{
					// actual cell = max(upper left cell+score, left cell-gap, upper cell-gap)
					double from_left(matrix.getValue(i,j-1)-gap);
					double from_above(matrix.getValue(i-1,j)-gap);
					double int1(spec1.getContainer()[i-1].getIntensity()), int2(spec2.getContainer()[j-1].getIntensity());
					double from_diagonal(matrix.getValue(i-1,j-1) + peakPairScore_(pos1,int1,pos2,int2,sigma));
					matrix.setValue(i,j,max(from_left,max(from_above,from_diagonal)));
				}
				else
				{
					// actual cell = max(left cell-gap, upper cell-gap)
					double from_left(matrix.getValue(i,j-1)-gap);
					double from_above(matrix.getValue(i-1,j)-gap);
					matrix.setValue(i,j,max(from_left,from_above));
				}
			}
		}
		
		#ifdef OPENMS_DEBUG 
		cout << endl << matrix << endl;
		#endif
		
		//get best overall score and return
		double best_score(DBL_MIN);
		for(UInt i = 0; i < matrix.cols(); i++)
		{
			best_score = max( best_score , matrix.getValue(matrix.rows()-1, i) );
		}
		for(UInt i = 0; i < matrix.rows(); i++)
		{
			best_score = max( best_score , matrix.getValue(i, matrix.cols()-1) );
		}
		
		//calculate selfalignment-scores for both input spectra
		double score_spec1(0),score_spec2(0);
		for (UInt i = 0; i < spec1.getContainer().size(); ++i)
		{
			double int_i(spec1.getContainer()[i].getIntensity());
			double pos_i(spec1.getContainer()[i].getMZ());
			score_spec1 += peakPairScore_(pos_i,int_i,pos_i,int_i,sigma);
		}
		for (UInt i = 0; i < spec2.getContainer().size(); ++i)
		{
			double int_i(spec2.getContainer()[i].getIntensity());
			double pos_i(spec2.getContainer()[i].getMZ());
			score_spec2 += peakPairScore_(pos_i,int_i,pos_i,int_i,sigma);
		}
		
		
		#ifdef OPENMS_DEBUG 
		cout << "score_spec1: " << score_spec1 << "score_spec2: " << score_spec2 << endl;
		#endif
		
		//normalize score to interval [0,1] with geometric mean
		double best_score_normalized(best_score/sqrt(score_spec1*score_spec2));
		
		/*
		#ifdef OPENMS_DEBUG 
		cout << "score_spec1: " << score_spec1 << " score_spec2: " << score_spec2 <<  " best_score: " << best_score << endl;
		#endif	
		
		//normalize score to interval [0,1] with arithmeic mean
		double best_score_normalized( (best_score*2) / (score_spec1 + score_spec2) );		
		*/
	    return best_score_normalized;
	
	}
	
	vector< pair<UInt,UInt> > PeakAlignment::getAlignmentTraceback (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const
	{
		const double epsilon = (double)param_.getValue("epsilon");
		
		//TODO gapcost dependence on distance ?
		const double gap = (double)param_.getValue("epsilon");

		//initialize alignment matrix with 0 in (0,0) and a multiple of gapcost in the first row/col matrix(row,col,values)
		Matrix<double> matrix(spec1.getContainer().size()+1, spec2.getContainer().size()+1,0);
		for (UInt i = 1; i < matrix.rows(); i++)
		{
			matrix.setValue(i,0,-gap*i);
		}
		for (UInt i = 1; i < matrix.cols(); i++)
		{
			matrix.setValue(0,i,-gap*i);
		}
		
		// gives the direction of the matrix cell that originated the respective cell
		// e.g. matrix(i+1,j+1) could have originated from matrix(i,j), matrix(i+1,j) or matrix(i,j+1)
		// so traceback(i,j) represents matrix(i+1,j+1) and contains a "1"-from diagonal, a "0"-from left or a "2"-from above
		Matrix<UInt> traceback	(spec1.getContainer().size(), spec2.getContainer().size());
		
		//get sigma - the standard deviation (sqrt of variance)
		double mid(0);
		for (UInt i = 0; i < spec1.getContainer().size(); ++i)
		{
			for (UInt j = 0; j < spec2.getContainer().size(); ++j)
			{
				double pos1(spec1.getContainer()[i].getMZ()), pos2(spec2.getContainer()[j].getMZ());
				mid += fabs(pos1 - pos2);
			}
		}
		mid = mid/(spec1.getContainer().size()*spec2.getContainer().size());
			
			#ifdef OPENMS_DEBUG	
			cout << mid << endl;
			#endif

		double var(0);
		for (UInt i = 0; i < spec1.getContainer().size(); ++i)
		{
			for (UInt j = 0; j < spec2.getContainer().size(); ++j)
			{
				double pos1(spec1.getContainer()[i].getMZ()), pos2(spec2.getContainer()[j].getMZ());
				var += (fabs(pos1 - pos2) - mid)*(fabs(pos1 - pos2) - mid);
			}
		}	
		var = var/(spec1.getContainer().size()*spec2.getContainer().size());
			
			#ifdef OPENMS_DEBUG	
			cout << var << endl;
			#endif

		const double sigma(sqrt(var));
	
			#ifdef OPENMS_DEBUG	
			cout << sigma << endl;
			#endif

		
		//fill alignment matrix
		for (UInt i = 1; i < spec1.getContainer().size()+1; ++i)
		{
			for (UInt j = 1; j < spec2.getContainer().size()+1; ++j)
			{
				double pos1(spec1.getContainer()[i-1].getMZ()), pos2(spec2.getContainer()[j-1].getMZ());
				//only if peaks are in reasonable proximity alignment is considered else only gaps
				if (fabs(pos1 - pos2) <= epsilon)
				{
					// actual cell = max(upper left cell+score, left cell-gap, upper cell-gap)
					double from_left(matrix.getValue(i,j-1)-gap);
					double from_above(matrix.getValue(i-1,j)-gap);
					double int1(spec1.getContainer()[i-1].getIntensity()), int2(spec2.getContainer()[j-1].getIntensity());
					double from_diagonal(matrix.getValue(i-1,j-1) + peakPairScore_(pos1,int1,pos2,int2,sigma));
					matrix.setValue(i,j,max(from_left,max(from_above,from_diagonal)));
					
					// TODO the cases where all or two values are equal
					if(from_diagonal>from_left && from_diagonal>from_above)
					{
						traceback.setValue(i-1,j-1,1);
					}
					else
					{
						if(from_left>from_diagonal && from_left>from_above)
						{
							traceback.setValue(i-1,j-1,0);
						}
						else
						{
							if(from_above>from_diagonal && from_above>from_left)
							{
								traceback.setValue(i-1,j-1,2);
							}
						}
					}
				}
				else
				{
					// actual cell = max(left cell-gap, upper cell-gap)
					double from_left(matrix.getValue(i,j-1)-gap);
					double from_above(matrix.getValue(i-1,j)-gap);
					matrix.setValue(i,j,max(from_left,from_above));
					if(from_left>from_above)
					{
						traceback.setValue(i-1,j-1,0);
					}
					else //from_left <= from_above 
					{
						traceback.setValue(i-1,j-1,2);
					}
				}
			}
		}
		//return track from best alloverscore to 0,0
		vector< pair<UInt,UInt> > ret_val;
		pair<UInt,UInt> max_pair();
		
		//get matrix coordinates from best alloverscore
		UInt row_index(0), col_index(0);
		double best_score(DBL_MIN);
		for(UInt i = 0; i < matrix.cols(); i++)
		{
			if(best_score < matrix.getValue(matrix.rows()-1, i) )
			{
				best_score = matrix.getValue(matrix.rows()-1, i);
				row_index = matrix.rows()-1;
				col_index = i;
			}
		}
		for(UInt i = 0; i < matrix.rows(); i++)
		{
			if(best_score < matrix.getValue(i, matrix.cols()-1) )
			{
				best_score = matrix.getValue(i, matrix.cols()-1) ;
				row_index = i;
				col_index = matrix.cols()-1;
			}
		}

		//TODO invariante prüfen
		while( row_index > 0 && col_index > 0 )
		{
			//from diagonal - peaks aligned
			if( traceback.getValue(row_index-1,col_index-1) == 1 )
			{
				//register aligned peaks only
				ret_val.insert(ret_val.begin(), pair<UInt,UInt>(row_index-1,col_index-1));
				row_index = row_index-1;
				col_index = col_index-1;
			}
			// gap alignment
			else if ( traceback.getValue(row_index-1,col_index-1) == 0 )
			{
				col_index = col_index-1;
			}
			else
			{
				row_index = row_index-1;
			}
		}		
		
		#ifdef OPENMS_DEBUG 
		cout << endl << matrix << endl << traceback << endl;
		#endif
		
	    return ret_val;	
	}

	
	double PeakAlignment::peakPairScore_(double& pos1, double& intens1, double& pos2, double& intens2, const double& sigma) const
	{
		//scoring formula : peakintensity score * peakposition score 
		double pi(sqrt(intens1*intens2));
		double pp( (1/(sigma*sqrt(2*M_PI)))* exp(-(fabs(pos1-pos2))/2*sigma*sigma) );
		
		#ifdef OPENMS_DEBUG 
		cout << fabs(pos1-pos2) << " - "<< pi*pp << endl;
		#endif
		
		return pi*pp;
	}
}
