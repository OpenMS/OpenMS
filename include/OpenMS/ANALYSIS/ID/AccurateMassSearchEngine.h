// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H
#define OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


#include <vector>

namespace OpenMS
{
    class OPENMS_DLLAPI AccurateMassSearchResult
    {
    public:
        /// Default constructor
        AccurateMassSearchResult();

        /// Default destructor
        virtual ~AccurateMassSearchResult();

        /// copy constructor
        AccurateMassSearchResult(const AccurateMassSearchResult& );

        /// assignment operator
        AccurateMassSearchResult & operator=(const AccurateMassSearchResult& );

        /// getter & setter methods
        DoubleReal getAdductMass()
        {
            return adduct_mass_;
        }

        DoubleReal getAdductMass() const
        {
            return adduct_mass_;
        }

        void setAdductMass(const DoubleReal& m)
        {
            adduct_mass_ = m;
        }

        DoubleReal getQueryMass()
        {
            return query_mass_;
        }

        void setQueryMass(const DoubleReal& m)
        {
            query_mass_ = m;
        }

        DoubleReal getFoundMass()
        {
            return found_mass_;
        }

        void setFoundMass(const DoubleReal& m)
        {
            found_mass_ = m;
        }

        DoubleReal getCharge()
        {
            return charge_;
        }

        DoubleReal getCharge() const
        {
            return charge_;
        }

        void setCharge(const DoubleReal& ch)
        {
            charge_ = ch;
        }

        DoubleReal getErrorPPM()
        {
            return error_ppm_;
        }

        void setErrorPPM(const DoubleReal& ppm)
        {
            error_ppm_ = ppm;
        }

        DoubleReal getObservedRT()
        {
            return observed_rt_;
        }

        DoubleReal getObservedRT() const
        {
            return observed_rt_;
        }

        void setObservedRT(const DoubleReal& rt)
        {
            observed_rt_ = rt;
        }

        DoubleReal getObservedIntensity()
        {
            return observed_intensity_;
        }

        DoubleReal getObservedIntensity() const
        {
            return observed_intensity_;
        }

        void setObservedIntensity(const DoubleReal& intensity)
        {
            observed_intensity_ = intensity;
        }

        DoubleReal getMatchingIndex()
        {
            return matching_index_;
        }

        void setMatchingIndex(const DoubleReal& idx)
        {
            matching_index_ = idx;
        }

        String getFoundAdduct()
        {
            return found_adduct_;
        }

        String getFoundAdduct() const
        {
            return found_adduct_;
        }

        void setFoundAdduct(const String& add)
        {
            found_adduct_ = add;
        }

        String getFormulaString()
        {
            return empirical_formula_;
        }

        String getFormulaString() const
        {
            return empirical_formula_;
        }

        void setEmpiricalFormula(const String& ep)
        {
            empirical_formula_ = ep;
        }

        std::vector<String> getMatchingHMDBids()
        {
            return matching_hmdb_ids_;
        }

        std::vector<String> getMatchingHMDBids() const
        {
            return matching_hmdb_ids_;
        }

        void setMatchingHMDBids(const std::vector<String>& match_ids)
        {
            matching_hmdb_ids_ = match_ids;
        }

        DoubleReal getIsotopesSimScore()
        {
            return isotopes_sim_score_;
        }

        DoubleReal getIsotopesSimScore() const
        {
            return isotopes_sim_score_;
        }

        void setIsotopesSimScore(const DoubleReal& sim_score)
        {
            isotopes_sim_score_ = sim_score;
        }

        // DoubleReal computeCombinedScore(); // not implemented
        // debug/output functions
        void outputResults() const;

    private:
        /// Stored information/results of DB query
        DoubleReal adduct_mass_;
        DoubleReal query_mass_;
        DoubleReal found_mass_;
        DoubleReal charge_;
        DoubleReal error_ppm_;
        DoubleReal observed_rt_;
        DoubleReal observed_intensity_;
        Size matching_index_;

        String found_adduct_;
        String empirical_formula_;
        std::vector<String> matching_hmdb_ids_;

        DoubleReal isotopes_sim_score_;
    };


    class OPENMS_DLLAPI AccurateMassSearchEngine :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    AccurateMassSearchEngine();

    // Explicit constructor
    AccurateMassSearchEngine(const String& map_fname);

    /// Default destructor
    virtual ~AccurateMassSearchEngine();

    void queryByMass(const DoubleReal&, const DoubleReal&, std::vector<AccurateMassSearchResult>&);
    void queryByFeature(const Feature&, std::vector<AccurateMassSearchResult>&);

    /// main method of AccurateMassSearchEngine
    void run(const FeatureMap<>&, MzTab&);


protected:
    virtual void updateMembers_();

private:
    /// private member functions
    void parseMappingFile_(const String&);
    void parseStructMappingFile_(const String&);
    void parseAdductsFile_(const String&);
    void searchMass_(const DoubleReal&, std::vector<Size>& hit_indices);

    void parseAdductString_(const String&, std::vector<String>&);
    void computeNeutralMassFromAdduct_(const DoubleReal&, const String&, DoubleReal&, DoubleReal&);

    DoubleReal computeCosineSim_(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y);
    DoubleReal computeEuclideanDist_(const std::vector<DoubleReal>& x, const std::vector<DoubleReal>& y);
    DoubleReal computeIsotopePatternSimilarity_(const Feature&, const EmpiricalFormula&);

    typedef std::vector<std::vector<AccurateMassSearchResult> > QueryResultsTable;

    void exportMzTab_(const QueryResultsTable&, MzTab&);

    /// private member variables
    typedef std::vector<std::vector<String> > MassIDMapping;
    typedef std::map<String, std::vector<String> > HMDBPropsMapping;

    std::vector<DoubleReal> masskey_table_;
    MassIDMapping mass_id_mapping_;
    HMDBPropsMapping hmdb_properties_mapping_;
    std::vector<String> mass_formula_mapping_;

    /// parameter stuff
    DoubleReal mass_error_value_;
    String mass_error_unit_;
    String ion_mode_;
    bool iso_similarity_;

    String pos_adducts_fname_;
    String neg_adducts_fname_;

    StringList pos_adducts_;
    StringList neg_adducts_;
  };


}




#endif // OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H
