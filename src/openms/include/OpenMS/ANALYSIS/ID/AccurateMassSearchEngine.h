// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H
#define OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/SYSTEM/File.h>

#include <vector>

namespace OpenMS
{
    class EmpiricalFormula;

    class OPENMS_DLLAPI AccurateMassSearchResult
    {
    public:
        /// Default constructor
        AccurateMassSearchResult();

        /// Default destructor
        ~AccurateMassSearchResult();

        /// copy constructor
        AccurateMassSearchResult(const AccurateMassSearchResult& );

        /// assignment operator
        AccurateMassSearchResult & operator=(const AccurateMassSearchResult& );

        /// getter & setter methods
        double getAdductMass() const;
        void setAdductMass(const double&);

        double getQueryMass() const;
        void setQueryMass(const double&);

        double getFoundMass() const;
        void setFoundMass(const double&);

        Int getCharge() const;
        void setCharge(const Int&);

        double getErrorPPM() const;
        void setErrorPPM(const double&);

        double getObservedRT() const;
        void setObservedRT(const double& rt);

        double getObservedIntensity() const;
        void setObservedIntensity(const double&);

        std::vector<double> getIndividualIntensities() const;
        void setIndividualIntensities(const std::vector<double>&);

        Size getMatchingIndex() const;
        void setMatchingIndex(const Size&);

        Size getSourceFeatureIndex() const;
        void setSourceFeatureIndex(const Size&);

        const String& getFoundAdduct() const;
        void setFoundAdduct(const String&);

        const String& getFormulaString() const;
        void setEmpiricalFormula(const String&);

        const std::vector<String>& getMatchingHMDBids() const;
        void setMatchingHMDBids(const std::vector<String>&);

        double getIsotopesSimScore() const;
        void setIsotopesSimScore(const double&);

        // double computeCombinedScore(); // not implemented
        // debug/output functions
        void outputResults() const;

    private:
        /// Stored information/results of DB query
        double adduct_mass_;
        double query_mass_;
        double found_mass_;
        Int charge_;
        double error_ppm_;
        double observed_rt_;
        double observed_intensity_;
        std::vector<double> individual_intensities_;
        Size matching_index_;
        Size source_feature_index_;

        String found_adduct_;
        String empirical_formula_;
        std::vector<String> matching_hmdb_ids_;

        double isotopes_sim_score_;
    };

  /**
    @brief An algorithm to search for exact mass matches from a spectrum against a database (e.g. HMDB).

    For each peak, neutral masses are reconstructed from observed (spectrum) m/z values by enumerating all possible adducts with matching charge.
    The resulting neutral masses (can be more than one, depending on list of possible adducts) are matched against masses from
    a database within a certain mass error (Da or ppm).

    Supports any database which contains an identifier, chemical sum formula and (optional) mass.
    If masses in the database are not given (= set to 0), they are computed from sum formulas.

    Both positive and negative ion mode is supported. Charge for (Consensus-)Features can be either positive or negative, but
    only the absolute value is used since many FeatureFinders will only report positive charges even in negative ion mode.
    Entities with charge=0 are treated as "unknown charge" and are tested with all potential adducts and subsequently matched against the database.

    A list of potential adducts can be given for each mode separately.

    Ionization mode of the observed m/z values can be determined automatically if the input map (either FeatureMap or ConsensusMap) is annotated
    with a meta value, as done by @ref TOPP_FeatureFinderMetabo.


    @ingroup Analysis_ID
  */
    class OPENMS_DLLAPI AccurateMassSearchEngine :
    public DefaultParamHandler,
    public ProgressLogger
  {
public:
    /// Default constructor
    AccurateMassSearchEngine();

    /// Default destructor
    virtual ~AccurateMassSearchEngine();

    /**
      @brief search for a specific observed mass by enumerating all possible adducts and search M+X against database

       */
    void queryByMass(const double& observed_mass, const Int& observed_charge, std::vector<AccurateMassSearchResult>& results);
    void queryByFeature(const Feature& feature, const Size& feature_index, std::vector<AccurateMassSearchResult>& results);
    void queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, std::vector<AccurateMassSearchResult>& results);

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    void run(FeatureMap&, MzTab&);

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    void run(ConsensusMap&, MzTab&);

    /// the internal ion-mode used depending on annotation of input data if "ion_mode" was set to 'auto'
    /// if run() was not called yet, this will be identical to 'ion_mode', i.e. 'auto' will no be resolved yet
    const String& getInternalIonMode();

protected:
    virtual void updateMembers_();

private:
    /// private member functions

    /// if ion-mode is auto, this will set the internal mode according to input data
    template <typename MAPTYPE> void resolveAutoMode_(const MAPTYPE& map)
    {
       String ion_mode_detect_msg = "";
       if (map.size() > 0 && map[0].metaValueExists("scan_polarity"))
       {
         StringList pols = ListUtils::create<String>(String(map[0].getMetaValue("scan_polarity")), ';');
         if (pols.size() == 1 && pols[0].size() > 0)
         {
            pols[0].toLower();
            if (pols[0] == "positive" || pols[0] == "negative")
            {
              ion_mode_internal_ = pols[0];
              LOG_INFO << "Setting auto ion-mode to '" << ion_mode_internal_ << "' for file " << File::basename(map.getLoadedFilePath()) << std::endl;
            }
            else ion_mode_detect_msg = String("Meta value 'scan_polarity' does not contain unknown ion mode") + String(map[0].getMetaValue("scan_polarity"));
         }
         else
         {
           ion_mode_detect_msg = String("ambiguous ion mode: ") + String(map[0].getMetaValue("scan_polarity"));
         }
       }
       else
       {
         ion_mode_detect_msg = String("Meta value 'scan_polarity' not found in first element of (Consensus-)Feature map!");
       }

       if (ion_mode_detect_msg.size() > 0)
       {
         throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("Auto ionization mode could not resolve ion mode of data (") + ion_mode_detect_msg + "!");
       }
    }

    /// parse database and adduct files
    void init_();

    bool is_initialized_; // true if init_() was called without any subsequent param changes

    void parseMappingFile_(const String&);
    void parseStructMappingFile_(const String&);
    void parseAdductsFile_(const String& filename, StringList& result);
    void searchMass_(const double&, std::vector<Size>& hit_indices);

    /// add search results to a Consensus/Feature
    void annotate_(const std::vector<AccurateMassSearchResult>&, BaseFeature&);

    /**
      @brief given an adduct and an observed mass, we compute the neutral mass (without adduct) and the theoretical charge (of the adduct)

    */
    void computeNeutralMassFromAdduct_(const double& observed_mass, const String& adduct_string, double& neutral_mass, Int& charge_value);

    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y);
    double computeEuclideanDist_(const std::vector<double>& x, const std::vector<double>& y);
    double computeIsotopePatternSimilarity_(const Feature&, const EmpiricalFormula&);


    String ion_mode_internal_;

    typedef std::vector<std::vector<AccurateMassSearchResult> > QueryResultsTable;

    void exportMzTab_(const QueryResultsTable&, MzTab&);

    /// private member variables
    typedef std::vector<std::vector<String> > MassIDMapping;
    typedef std::map<String, std::vector<String> > HMDBPropsMapping;

    struct MappingEntry_
    {
      double mass;
      std::vector<String> massIDs;
      String formula;
    };
    std::vector<MappingEntry_> mass_mappings_;

    struct CompareEntryAndMass_     // defined here to allow for inlining by compiler
    {
       double asMass( const MappingEntry_& v ) const
       {
          return v.mass;
       }

       double asMass( double t ) const
       {
          return t;
       }

       template< typename T1, typename T2 >
       bool operator()( T1 const& t1, T2 const& t2 ) const
       {
           return asMass(t1) < asMass(t2);
       }
    };

    HMDBPropsMapping hmdb_properties_mapping_;

    /// parameter stuff
    double mass_error_value_;
    String mass_error_unit_;
    String ion_mode_;
    bool iso_similarity_;

    String pos_adducts_fname_;
    String neg_adducts_fname_;

    String db_mapping_file_;
    String db_struct_file_;

    StringList pos_adducts_;
    StringList neg_adducts_;
  };

}

#endif // OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H
