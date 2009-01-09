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
//

#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <math.h>
#include <string>
	
using namespace OpenMS;
using namespace std;


	PeakIntensityPredictor::PeakIntensityPredictor()
		: llm_()
	{
	}

	PeakIntensityPredictor::~PeakIntensityPredictor()
	{
	}

	DoubleReal PeakIntensityPredictor::predict(const AASequence& sequence)
	{
		//get corresponding attribute values
		vector<DoubleReal> aafeat = AAIndex::getPropertyVector(sequence);
		//normalize vector (center and scale by variance)
		llm_.normalizeVector(aafeat);
		//pass data_ to llm model and return predicted value
		return map_(aafeat);
	}

	DoubleReal PeakIntensityPredictor::predict(const AASequence& sequence, vector<DoubleReal>& add_info)
	{
		//get corresponding attribute values
		vector<DoubleReal> aafeat = AAIndex::getPropertyVector(sequence);
		//normalize vector (center and scale by variance)
		llm_.normalizeVector(aafeat);
		//pass data_ to llm model and return predicted value
		DoubleReal tmp = map_(aafeat);

		//Determine corresponding winning prototype and number of cluster
		//peptide is assigned to. Additionally, the error is given.	
		add_info = calculateAddInfo_(aafeat);
		
		return tmp;
	}

	vector<DoubleReal> PeakIntensityPredictor::predict(const vector<AASequence>& sequences)
	{
		vector<DoubleReal> out(sequences.size());

		for(Size i=0; i<sequences.size(); i++)
		{
			out[i] = predict(sequences[i]);
		}
		
		return out;
	}

	vector<DoubleReal> PeakIntensityPredictor::predict(const vector<AASequence>& sequences, vector<vector<DoubleReal> >& add_info)
	{
		vector<DoubleReal> out(sequences.size());
		add_info.resize(sequences.size());

		//for each element in sequences
		for(Size i=0; i<sequences.size(); i++)
		{
			out[i] = predict(sequences[i],add_info[i]);
		}
		
		return out;
	}

	DoubleReal PeakIntensityPredictor::map_(const vector<DoubleReal>& data)
	{
		DoubleReal c_x = 0.0;
		DoubleReal sum_g_i, g_i;
		Matrix<DoubleReal> code = llm_.getCodebooks();
		vector<DoubleReal> wout = llm_.getVectorWout();
		Matrix<DoubleReal> A = llm_.getMatrixA();

		//determine best matching unit
		UInt winner = findWinner_(data);
		//calculate Gaussian neighborhood function 
		vector< DoubleReal > nei = llm_.neigh(llm_.getCord(), winner, llm_.getLLMParam().radius);

		sum_g_i = 0.0;
		//for each prototype 
		for(Size c=0; c<code.rows(); c++)
		{
			//sum up total sum of nei
			sum_g_i += nei[c];
		}
		g_i = 0.0;
			
		//map data to best matching unit of prototypes and get 
		//linear mapping 
		for (Size i = 0; i < code.rows(); i++) 
		{
			c_x = 0.0;
			for(Size c=0; c<code.cols(); c++)
			{
				c_x += (A.getValue(i,c) * (data[c] - code.getValue(i, c)));
			}
			//add linear bias to mapped value
			g_i += (c_x + wout[i])*nei[i];
		}
		//divide by total neighborhod values
		g_i = g_i/sum_g_i;
		
		//normalize predicted values to distribution of training data
		g_i -= 3.364288; //mean(targets)
		g_i /= 1.332298; //sd(targets)
		
		return g_i;
	}

	UInt PeakIntensityPredictor::findWinner_(const vector<DoubleReal>& data)
	{
		UInt winner = 0;
		DoubleReal min_dist = 0.0;
		Matrix<DoubleReal> code = llm_.getCodebooks();

		//calculate euclidean distance of vector data to prototype no 1.
		for(Size c=0; c<data.size(); c++)
		{
			min_dist += (data[c] - code.getValue(1, c))*(data[c] - code.getValue(1, c));
		}
	
		//calculate euclidean distance of vector data to the remaining prototypes
		for (Size i = 1; i < code.rows(); i++) 
		{
			DoubleReal dd = 0.0;
			for(Size c=0; c<data.size(); c++)
			{
				dd += (data[c] - code.getValue(i, c))*(data[c] - code.getValue(i, c));
			}
			//which is min?
			if (dd < min_dist) 
			{
				winner = i;
				min_dist = dd;
			}
		}
		
		//return number of winning prototype with minimum distance to data
		return winner;
	}

	vector<DoubleReal> PeakIntensityPredictor::calculateAddInfo_(const vector<DoubleReal>& data) 
	{
		vector<DoubleReal> foo(3);
		UInt winner = findWinner_(data);
		Matrix<DoubleReal> code = llm_.getCodebooks();
		Matrix<UInt> cord = llm_.getCord();

		foo[0] = cord.getValue(winner,0);
		foo[1] = cord.getValue(winner,1);

		DoubleReal dd = 0.0;
		//get distance of data to best matching unit (winner).
		for(Size c = 0; c < data.size(); c++)
		{
			dd += (data[c] - code.getValue(winner,c))*(data[c] - code.getValue(winner,c));
		}
		//store the minimum distance in third column 
		foo[2] = sqrt(dd);
		
		return foo;
	}










