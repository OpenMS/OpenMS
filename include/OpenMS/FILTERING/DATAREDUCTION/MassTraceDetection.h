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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#ifndef OPENMS_MASSTRACEDETECTION_H
#define OPENMS_MASSTRACEDETECTION_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
    class OPENMS_DLLAPI MassTraceDetection :
            public DefaultParamHandler,
            public ProgressLogger
    {
    public:
        /// Default constructor
        MassTraceDetection();

        /// Default destructor
        virtual ~MassTraceDetection();


        /// main method of MassTraceDetection
        void run(const MSExperiment<Peak1D>&, std::vector<MassTrace>&);





    protected:
        virtual void updateMembers_();


    private:
        /// parameter stuff
        DoubleReal mass_error_ppm_;
        DoubleReal noise_threshold_int_;
        DoubleReal chrom_apex_snt_;
        DoubleReal chrom_fwhm_;

        DoubleReal min_sample_rate_;

    };


}




#endif // OPENMS_MASSTRACEDETECTION_H
