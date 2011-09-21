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
// $Maintainer: $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#ifndef LowessSmoothing_H
#define LowessSmoothing_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
// #include <OpenMS/CONCEPT/ProgressLogger.h>


namespace OpenMS
{
    class OPENMS_DLLAPI LowessSmoothing
        : public DefaultParamHandler
    {
    public:
        /// Constructor
        LowessSmoothing();

        /// Destructor
        virtual ~LowessSmoothing();

        /**

            @brief LOWESS (locally weighted scatterplot smoothing) filtering. Fitting simple models to localized subsets of the data to build up a function that describes the deterministic part of the variation in the data, point by point.

        */
        typedef std::vector<DoubleReal> DoubleVector;

        void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&);



    protected:
        virtual void updateMembers_();


    private:
        DoubleReal window_size_;

        DoubleReal tricube_(DoubleReal, DoubleReal);

    };


} // namespace OpenMS
#endif // LowessSmoothing_H
