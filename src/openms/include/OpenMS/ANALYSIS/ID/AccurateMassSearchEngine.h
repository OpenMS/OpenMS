// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
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
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

#include <iosfwd>
#include <vector>

namespace OpenMS
{ 
  class EmpiricalFormula;

  class OPENMS_DLLAPI AdductInfo
  {

  public: 
    /**
      C'tor, to build a representation of an adduct.

      @param name Identifier as given in the Positive/Negative-Adducts file, e.g. 'M+2K-H;1+'
      @param adduct Formula of the adduct, e.g. '2K-H'
      @param charge The charge (must not be 0; can be negative), e.g. 1
      @param is_intrinsic True for a molecule without an explicit adduct, e.g. 'M;-1'
      @param mol_multiplier Molecular multiplier, e.g. for charged dimers '2M+H;+1'

    **/
    AdductInfo(const String& name, const EmpiricalFormula& adduct, int charge, UInt mol_multiplier = 1);

    /// returns the neutral mass of the small molecule without adduct (creates monomer from nmer, decharges and removes the adduct (given m/z of [nM+Adduct]/|charge| returns mass of [M])
    double getNeutralMass(double observed_mz) const;

    /// returns the m/z of the small molecule with neutral mass @p neutral_mass if the adduct is added (given mass of [M] returns m/z of [nM+Adduct]/|charge|)
    double getMZ(double neutral_mass) const;

    /// checks if an adduct (e.g.a 'M+2K-H;1+') is valid, i.e if the losses (==negative amounts) can actually be lost by the compound given in @p db_entry.
    /// If the negative parts are present in @p db_entry, true is returned.
    bool isCompatible(EmpiricalFormula db_entry) const;

    /// get charge of adduct
    int getCharge() const;

    /// original string used for parsing
    const String& getName() const;

    /// parse an adduct string containing a formula (must contain 'M') and charge, separated by ';'.
    /// e.g. M+H;1+
    /// 'M' can have multipliers, e.g. '2M + H;1+' (for a singly charged dimer)
    static AdductInfo parseAdductString(const String& adduct);

  private:
    /// hide default C'tor
    AdductInfo();

    /// members
    String name_; ///< arbitrary name, only used for error reporting
    EmpiricalFormula ef_; ///< EF for the actual adduct e.g. 'H' in 2M+H;+1
    double mass_; ///< computed from ef_.getMonoWeight(), but stored explicitly for efficiency
    int charge_;  ///< negative or positive charge; must not be 0
    UInt mol_multiplier_; ///< Mol multiplier, e.g. 2 in 2M+H;+1
  };

  class OPENMS_DLLAPI AccurateMassSearchResult
  {
  public:
    /// Default constructor
    AccurateMassSearchResult();

    /// Default destructor
    ~AccurateMassSearchResult();

    /// copy constructor
    AccurateMassSearchResult(const AccurateMassSearchResult&);

    /// assignment operator
    AccurateMassSearchResult& operator=(const AccurateMassSearchResult&);

    /// get the m/z of the small molecule + adduct
    double getObservedMZ() const;

    /// set the m/z of the small molecule + adduct
    void setObservedMZ(const double&);

    /// get the theoretical m/z of the small molecule + adduct
    double getCalculatedMZ() const;

    /// set the theoretical m/z of the small molecule + adduct
    void setCalculatedMZ(const double&);

    /// get the mass used to query the database (uncharged small molecule)
    double getQueryMass() const;

    /// set the mass used to query the database (uncharged small molecule)
    void setQueryMass(const double&);

    /// get the mass returned by the query (uncharged small molecule)
    double getFoundMass() const;

    /// set the mass returned by the query (uncharged small molecule)
    void setFoundMass(const double&);

    /// get the charge
    Int getCharge() const;

    /// set the charge
    void setCharge(const Int&);

    /// get the error between observed and theoretical m/z in ppm
    double getMZErrorPPM() const;

    /// set the error between observed and theoretical m/z in ppm
    void setMZErrorPPM(const double);

    /// get the observed rt
    double getObservedRT() const;

    /// set the observed rt
    void setObservedRT(const double& rt);

    /// get the observed intensity
    double getObservedIntensity() const;

    /// set the observed intensity
    void setObservedIntensity(const double&);

    /// get the observed intensities
    std::vector<double> getIndividualIntensities() const;

    /// set the observed intensities
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

    /// return trace intensities of the underlying feature;
    const std::vector<double>& getMasstraceIntensities() const;
    void setMasstraceIntensities(const std::vector<double>&);
    
    double getIsotopesSimScore() const;
    void setIsotopesSimScore(const double&);

    // debug/output functions
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr);

private:
    /// Stored information/results of DB query
    double observed_mz_;
    double theoretical_mz_;
    double searched_mass_;
    double db_mass_;
    Int charge_;
    double mz_error_ppm_;
    double observed_rt_;
    double observed_intensity_;
    std::vector<double> individual_intensities_;
    Size matching_index_;
    Size source_feature_index_;

    String found_adduct_;
    String empirical_formula_;
    std::vector<String> matching_hmdb_ids_;
    
    std::vector<double> mass_trace_intensities_;
    double isotopes_sim_score_;
  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const AccurateMassSearchResult& amsr);

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

    A file with a list of potential adducts can be given for each mode separately. 
    Each line contains a chemical formula (plus quantor) and a charge (separated by semicolon), e.g.
    M+H;1+ 
    The M can be preceded by a quantor (e.g.2M, 3M), implicitly assumed as 1.
    The chemical formula can contain multiple segments, separated by + or - operators, e.g. M+H-H2O;+1 (water loss in positive mode).
    Brackets are implicit per segment, i.e. M+H-H2O is parsed as M + (H) - (H2O).
    Each segment can also be preceded by a quantor, e.g. M+H-H2O would parse as
    M + (H) - 2x(H2O).
    If debug mode is enabled, the masses of each segment are printed for verification.
    In particular, typing H20 (twenty H) is different from H2O (water).

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
    ~AccurateMassSearchEngine() override;

    /**
      @brief search for a specific observed mass by enumerating all possible adducts and search M+X against database

       */
    void queryByMZ(const double& observed_mz, const Int& observed_charge, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const;
    void queryByFeature(const Feature& feature, const Size& feature_index, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const;
    void queryByConsensusFeature(const ConsensusFeature& cfeat, const Size& cf_index, const Size& number_of_maps, const String& ion_mode, std::vector<AccurateMassSearchResult>& results) const;

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    void run(FeatureMap&, MzTab&) const;

    /// main method of AccurateMassSearchEngine
    /// input map is not const, since it will get annotated with results
    /// @note Call init() before calling run!
    void run(ConsensusMap&, MzTab&) const;

    /// parse database and adduct files
    void init();

protected:
    void updateMembers_() override;

private:
    /// private member functions

    /// if ion-mode is auto, this will set the internal mode according to input data
    /// @throw InvalidParameter if ion mode cannot be resolved
    template <typename MAPTYPE> String resolveAutoMode_(const MAPTYPE& map) const
    {
      String ion_mode_internal;
      String ion_mode_detect_msg = "";
      if (map.size() > 0)
      {
        if (map[0].metaValueExists("scan_polarity"))
        {
          StringList pols = ListUtils::create<String>(String(map[0].getMetaValue("scan_polarity")), ';');
          if (pols.size() == 1 && pols[0].size() > 0)
          {
            pols[0].toLower();
            if (pols[0] == "positive" || pols[0] == "negative")
            {
              ion_mode_internal = pols[0];
              LOG_INFO << "Setting auto ion-mode to '" << ion_mode_internal << "' for file " << File::basename(map.getLoadedFilePath()) << std::endl;
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
          ion_mode_detect_msg = String("Meta value 'scan_polarity' not found in (Consensus-)Feature map");
        }
      }
      else
      { // do nothing, since map is
        LOG_INFO << "Meta value 'scan_polarity' cannot be determined since (Consensus-)Feature map is empty!" << std::endl;
      }

      if (ion_mode_detect_msg.size() > 0)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Auto ionization mode could not resolve ion mode of data (") + ion_mode_detect_msg + "!");
      }

      return ion_mode_internal;
    }
    
    void parseMappingFile_(const StringList&);
    void parseStructMappingFile_(const StringList&);
    void parseAdductsFile_(const String& filename, std::vector<AdductInfo>& result);
    void searchMass_(double neutral_query_mass, double diff_mass, std::pair<Size, Size>& hit_indices) const;

    /// add search results to a Consensus/Feature
    void annotate_(const std::vector<AccurateMassSearchResult>&, BaseFeature&) const;

    /// For two vectors of identical length, compute the cosine of the angle between them.
    /// Since we look at the angle, scaling of the vectors does not change the result (when ignoring numerical instability).
    double computeCosineSim_(const std::vector<double>& x, const std::vector<double>& y) const;

    double computeIsotopePatternSimilarity_(const Feature& feat, const EmpiricalFormula& form) const;

    typedef std::vector<std::vector<AccurateMassSearchResult> > QueryResultsTable;

    void exportMzTab_(const QueryResultsTable& overall_results, const Size number_of_maps, MzTab& mztab_out) const;

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

    struct CompareEntryAndMass_ // defined here to allow for inlining by compiler
    {
      double asMass(const MappingEntry_& v) const
      {
        return v.mass;
      }

      double asMass(double t) const
      {
        return t;
      }

      template <typename T1, typename T2>
      bool operator()(T1 const& t1, T2 const& t2) const
      {
        return asMass(t1) < asMass(t2);
      }

    };

    HMDBPropsMapping hmdb_properties_mapping_;

    bool is_initialized_; ///< true if init_() was called without any subsequent param changes

    /// parameter stuff
    double mass_error_value_;
    String mass_error_unit_;
    String ion_mode_;
    bool iso_similarity_;

    String pos_adducts_fname_;
    String neg_adducts_fname_;

    StringList db_mapping_file_;
    StringList db_struct_file_;

    std::vector<AdductInfo> pos_adducts_;
    std::vector<AdductInfo> neg_adducts_;

    String database_name_;
    String database_version_;

    bool keep_unidentified_masses_;
  };

}

#endif // OPENMS_ANALYSIS_ID_ACCURATEMASSSEARCHENGINE_H
