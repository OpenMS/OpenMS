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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
	namespace Math
	{
    Real pearsonCorrelation(std::vector<Real> model, std::vector<Real> data)
		{
			if (model.size()!=data.size())
			{
				//TODO throw exception
			}
			
			//calculate averages
			Real model_avg;
			std::accumulate(model.begin(),model.end(),model_avg);
			model_avg /= model.size();
			Real data_avg;
			std::accumulate(data.begin(),data.end(),data_avg);
			data_avg /= data.size();

			Real cross_product_sum = 0;
			Real data_square_sum   = 0;
			Real model_square_sum  = 0;
			for (UInt i=0; i<data.size();++i)
			{
				cross_product_sum += ( model[i] - model_avg) * ( data[i] - data_avg);
				data_square_sum += ( data[i] - data_avg)  * ( data[i] - data_avg);
				model_square_sum += ( model[i] - model_avg)  * ( model[i] - model_avg);			
			}

			if ( ! data_square_sum || ! model_square_sum ) return 0;
			
			return cross_product_sum / sqrt(data_square_sum * model_square_sum);
		}

	} // namespace Math
} // namespace OpenMS

