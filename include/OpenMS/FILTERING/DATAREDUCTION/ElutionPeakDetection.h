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


#ifndef OPENMS_ElutionPeakDetection_H
#define OPENMS_ElutionPeakDetection_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MassTrace.h>


namespace OpenMS
{

    class OPENMS_DLLAPI ElutionPeakDetection
        : public DefaultParamHandler, public ProgressLogger
    {
    public:
        /// Constructor
        ElutionPeakDetection();

        /// Destructor
        virtual ~ElutionPeakDetection();

        void detectPeaks(MassTrace&, std::vector<MassTrace>&);
        void detectPeaks(std::vector<MassTrace>&, std::vector<MassTrace>&);

        void filterByPeakWidth(std::vector<MassTrace>&, std::vector<MassTrace>&, DoubleReal);


    protected:
        virtual void updateMembers_();

    private:
        // DoubleReal chrom_fwhm_;
        Size window_size_;
        DoubleReal sample_rate_;

        String pw_filtering_;

        void detectElutionPeaks_(MassTrace&, std::vector<MassTrace>&);



    };
} // namespace OpenMS
#endif // OPENMS_ElutionPeakDetection_H
