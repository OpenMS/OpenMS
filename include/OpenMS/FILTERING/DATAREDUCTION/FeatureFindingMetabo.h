// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <svm.h>

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

        if (iso_pattern_.size() > 0)
        {
            label = iso_pattern_[0]->getLabel();
        }

        for (Size i = 1; i < iso_pattern_.size(); ++i)
        {
            String tmp_str = "_" + iso_pattern_[i]->getLabel();
            label += tmp_str;
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

    std::vector<DoubleReal> getAllIntensities(bool smoothed = false)
    {
        std::vector<DoubleReal> tmp;

        for (Size i = 0; i < iso_pattern_.size(); ++i)
        {
            if (!smoothed)
            {
                tmp.push_back(iso_pattern_[i]->getIntensity(false));
            }
            else
            {
                tmp.push_back(iso_pattern_[i]->getIntensity(true));
            }

        }

        return tmp;
    }


    DoubleReal getCentroidMZ()
    {
        if (iso_pattern_.empty())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no centroid MZ!", String(iso_pattern_.size()));
        }

        return iso_pattern_[0]->getCentroidMZ();
    }

    DoubleReal getCentroidRT()
    {
        if (iso_pattern_.empty())
        {
            throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "FeatureHypothesis is empty, no centroid RT!", String(iso_pattern_.size()));
        }

        iso_pattern_[0]->updateWeightedMeanRT();

        return iso_pattern_[0]->getCentroidRT();
    }

    DoubleReal getFWHM(bool use_smoothed_ints = false)
    {
        if (iso_pattern_.empty())
        {
            return 0.0;
        }

        return iso_pattern_[0]->estimateFWHM(use_smoothed_ints);
    }


    /// addMassTrace
    void addMassTrace(MassTrace&);
    DoubleReal getMonoisotopicFeatureIntensity(bool);
    DoubleReal getSummedFeatureIntensity(bool);


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
    DoubleReal computeOLSCoeff(const std::vector<DoubleReal>&, const std::vector<DoubleReal>&);
    DoubleReal computeCosineSim(const std::vector<DoubleReal>&, const std::vector<DoubleReal>&);

    svm_model *isotope_filt_svm;
    std::vector<DoubleReal> svm_feat_centers;
    std::vector<DoubleReal> svm_feat_scales;
    bool isLegalIsotopePattern_(FeatureHypothesis&);
    //bool isLegalAveraginePattern(FeatureHypothesis&);
    void loadIsotopeModel_();

    DoubleReal scoreMZ_(const MassTrace&, const MassTrace&, Size, Size);
    DoubleReal scoreRT_(const MassTrace&, const MassTrace&);

    DoubleReal computeAveragineSimScore(const std::vector<DoubleReal>&, const DoubleReal&);

    // DoubleReal scoreTraceSim_(MassTrace, MassTrace);
    // DoubleReal scoreIntRatio_(DoubleReal, DoubleReal, Size);
    void findLocalFeatures_(std::vector<MassTrace*>&, std::vector<FeatureHypothesis>&);


    /// parameter stuff
    DoubleReal local_rt_range_;
    DoubleReal local_mz_range_;
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    //DoubleReal mass_error_ppm_;
    DoubleReal chrom_fwhm_;

    bool report_summed_ints_;
    bool disable_isotope_filtering_;
    String isotope_model_;
    bool use_smoothed_intensities_;

};


}




#endif // OPENMS_FILTERING_DATAREDUCTION_FEATUREFINDINGMETABO_H
