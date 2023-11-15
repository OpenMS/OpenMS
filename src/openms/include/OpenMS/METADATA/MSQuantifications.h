// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/METADATA/MetaInfo.h>

#include <vector>
#include <map>

namespace OpenMS
{
  class OPENMS_DLLAPI MSQuantifications :
    public ExperimentalSettings
  {
public:
    /// @name Base type definitions
    //@{
    /// typedef docu
    typedef CVTermList ParamGroupList;         // userparams are exclusively inside the CVTermList's MetaInfoInterface

    enum QUANT_TYPES {MS1LABEL = 0, MS2LABEL, LABELFREE, SIZE_OF_QUANT_TYPES};       // derived from processing applied
    static const std::string NamesOfQuantTypes[SIZE_OF_QUANT_TYPES];

    //@}
    //~ InputFiles: //~ searchdb abbildung version,releasedate,#entries,dbname ber paramgrouplist
    //~ struct ParamGroupList
    //~ {
    //~ ParamGroupList()
    //~ {
    //~ }

    //~ ParamGroupList(const ParamGroupList& rhs)
    //~ :	cv_params(rhs.cv_params)
    //~ {
    //~ }

    //~ ~ParamGroupList()
    //~ {
    //~ }

    //~ ParamGroupList& operator = (const ParamGroupList& rhs)
    //~ {
    //~ if (&rhs != this)
    //~ {
    //~ cv_params = rhs.cv_params;
    //~ user_params = rhs.user_params;
    //~ }
    //~ return *this;
    //~ }

    //~ MetaInfoInterface user_params;
    //~ CVTermList cv_params;
    //~ };

    struct AnalysisSummary
    {
      AnalysisSummary() = default;

      AnalysisSummary(const AnalysisSummary& rhs) = default;

      virtual ~AnalysisSummary() = default;

      AnalysisSummary& operator=(const AnalysisSummary& rhs) = default;

      MetaInfo user_params_;
      CVTermList cv_params_;
      QUANT_TYPES quant_type_ = SIZE_OF_QUANT_TYPES;
    };

    struct Assay
    {
      //TODO feature_maps_ also in Assay?! srsly?!
      Assay() = default;

      Assay(const Assay& rhs) = default;

      virtual ~Assay() = default;

      Assay& operator=(const Assay& rhs) = default;

      String uid_;
      std::vector<std::pair<String, double> > mods_;
      std::vector<ExperimentalSettings> raw_files_;
      std::map<size_t, FeatureMap > feature_maps_;           // iTRAQ needs no FeatureMaps so ExperimentalSettings are not directly mapped to FeatureMaps
    };

    // TODO handle referencing from consensusmaps to featuremaps/rawfiles
    // TODO add ContactPerson or something to (Consensus)FeatureMap or DataProcessing (see below)
    // TODO rewrite OpenMS::DataProcessing - data not yet linked in openms core formats - below should go in analysissummary of MSQuantifications - input/output not possible to be carried along
    //~ if(DataProcessing::NamesOfProcessingAction[*it] == String("Quantitation"))
    //~ {
    //~ if (processing.getSoftware().getName()==String("SILACAnalyzer"))
    //~ {
    //~ experiment_type = MS1LABEL;
    //~ }
    //~ else if (processing.getSoftware().getName()==String("ITRAQAnalyzer"))
    //~ {
    //~ experiment_type = MS2LABEL;
    //~ }
    //~ else
    //~ {
    //~ experiment_type = LABELFREE;
    //~ }
    //~ }
    //~ QUANT_TYPES experiment_type = MS1LABEL;

    /// Constructor
    MSQuantifications() = default;

    /// Detailed Constructor
    MSQuantifications(const FeatureMap& fm, 
                      ExperimentalSettings& es,
                      std::vector<DataProcessing>& dps, std::vector<std::vector<std::pair<String, double>>> labels = {});

    /// Copy constructor
    MSQuantifications(const MSQuantifications & source) = default;

    /// Move constructor
    MSQuantifications(MSQuantifications&&) = default;

    /// Destructor
    ~MSQuantifications() override;

    /// Assignment operator
    MSQuantifications & operator=(const MSQuantifications & source) = default;

    /// Move assignment operator
    MSQuantifications& operator=(MSQuantifications&&) & = default;

    /// Equality operator
    bool operator==(const MSQuantifications & rhs) const;

    /// Equality operator
    bool operator!=(const MSQuantifications & rhs) const;

    /**
        @brief Loads data from a text file.

        @param filename The input file name.
        @param trim_lines Whether or not the lines are trimmed when reading them from file.
        @param first_n If set, only @p first_n lines the lines from the beginning of the file are read.

        @note this function uses unix-style line breaks

        @exception Exception::FileNotFound is thrown if the file could not be opened.

        TODO : implement
    */
    // void load(const String & filename, bool trim_lines = false, Int first_n = -1);

    const std::vector<DataProcessing> getDataProcessingList() const;
    const std::vector<Assay> & getAssays() const;
    std::vector<Assay> & getAssays();
    // std::map<String, ConsensusFeature::Ratio> & getRatios(); // TODO : implement
    const std::vector<ConsensusMap> & getConsensusMaps() const;
    std::vector<ConsensusMap> & getConsensusMaps();
    void setConsensusMaps(const std::vector<ConsensusMap> & );
    const std::vector<FeatureMap > & getFeatureMaps() const;
    const AnalysisSummary & getAnalysisSummary() const;
    AnalysisSummary & getAnalysisSummary();
    void setDataProcessingList(const std::vector<DataProcessing> & dpl);
    void setAnalysisSummaryQuantType(QUANT_TYPES r);
    void addConsensusMap(ConsensusMap & m);
    void assignUIDs();
    void registerExperiment(PeakMap & exp, std::vector<std::vector<std::pair<String, double> > > labels);
    void registerExperiment(ExperimentalSettings & es, std::vector<DataProcessing>& dp, std::vector<std::vector<std::pair<String, double> > > labels = (std::vector<std::vector<std::pair<String, double> > >()));

private:
    AnalysisSummary analysis_summary_;
    std::vector<MetaInfo> bibliographic_reference_;
    std::vector<ConsensusMap> consensus_maps_;
    std::vector<FeatureMap > feature_maps_;
    std::vector<Assay> assays_;
    std::vector<DataProcessing> data_processings_;
    //~ std::map<String,ConsensusFeature::Ratio > ratio_calculations_;
  };

} // namespace OpenMS

