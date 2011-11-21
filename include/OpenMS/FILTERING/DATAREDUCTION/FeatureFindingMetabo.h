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

#ifndef OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
#define OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <vector>

namespace OpenMS
{
    class OPENMS_DLLAPI CmpMassTraceByMZ
    {
    public:

        bool operator()(MassTrace x, MassTrace y) const
        {
            return x.getCentroidMZ() < y.getCentroidMZ();
        }
    };


    class OPENMS_DLLAPI FeatureHypothesis
    {
    public:
        /// default constructor
        FeatureHypothesis();

        /// default destructor
        ~FeatureHypothesis();

        /// copy constructor
        FeatureHypothesis(const FeatureHypothesis&);

        /// assignment operator
        FeatureHypothesis& operator=(const FeatureHypothesis& rhs);


        // getter & setter
        Size getSize() const
        {
            return iso_pattern_.size();
        }


        String getLabel()
        {
            String label;

            for (Size i = 0; i < iso_pattern_.size(); ++i)
            {
                label += iso_pattern_[i]->getLabel();
            }

            return label;
        }

        std::vector<String> getLabels()
        {
            std::vector<String> tmp_labels;

            for (Size i = 0; i < iso_pattern_.size(); ++i)
            {
                tmp_labels.push_back(iso_pattern_[i]->getLabel());
            }

            return tmp_labels;
        }

        DoubleReal getScore()
        {
            return feat_score_;
        }

        void setScore(const DoubleReal& score)
        {
            feat_score_ = score;
        }

        SignedSize getCharge()
        {
            return charge_;
        }

        void setCharge(const SignedSize& ch)
        {
            charge_ = ch;
        }


        DoubleReal getCentroidMZ()
        {
            if (iso_pattern_.size() == 0)
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no centroid MZ!", String(iso_pattern_.size()));
            }

            return iso_pattern_[0]->getCentroidMZ();
        }

        DoubleReal getCentroidRT()
        {
            if (iso_pattern_.size() == 0)
            {
                throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no centroid RT!", String(iso_pattern_.size()));
            }

            return iso_pattern_[0]->getSmoothedMaxRT();
        }

        DoubleReal getFWHM()
        {
            if (iso_pattern_.size() == 0)
            {
                return 0.0;
            }

            return iso_pattern_[0]->estimateFWHM();
        }


        /// addMassTrace
        void addMassTrace(MassTrace&);
        DoubleReal computeFeatureIntensity();


        Size getNumFeatPoints() const;
        std::vector<ConvexHull2D> getConvexHulls() const;

    private:
        // pointers of MassTraces contained in isotopic pattern
        std::vector<MassTrace*> iso_pattern_;
        DoubleReal feat_score_;

        SignedSize charge_;

    };


    class OPENMS_DLLAPI CmpHypothesesByScore
    {
    public:

        bool operator()(FeatureHypothesis x, FeatureHypothesis y) const
        {
            return x.getScore() > y.getScore();
        }
    };



    class OPENMS_DLLAPI FeatureFindingMetabo :
            public DefaultParamHandler,
            public ProgressLogger
    {
    public:
        /// Default constructor
        FeatureFindingMetabo();

        /// Default destructor
        virtual ~FeatureFindingMetabo();


        /// main method of FeatureFindingMetabo
        void run(std::vector<MassTrace>&, FeatureMap<>&);





    protected:
        virtual void updateMembers_();


    private:
        /// private member functions
        DoubleReal scoreMZ_(DoubleReal, DoubleReal, Size, Size);
        DoubleReal scoreRT_(DoubleReal, DoubleReal);
        DoubleReal scoreTraceSim_(MassTrace, MassTrace);
        DoubleReal scoreIntRatio_(DoubleReal, DoubleReal, Size);
        void findLocalFeatures_(std::vector<MassTrace*>&, std::vector<FeatureHypothesis>&);


        /// parameter stuff
        DoubleReal local_rt_range_;
        DoubleReal local_mz_range_;
        Size charge_lower_bound_;
        Size charge_upper_bound_;
        DoubleReal mass_error_ppm_;
        DoubleReal chrom_fwhm_;

    };


}




#endif // OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
