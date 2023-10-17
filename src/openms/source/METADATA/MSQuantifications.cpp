// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MSQuantifications.h>

#include <OpenMS/CONCEPT/UniqueIdGenerator.h>

#include <utility>

using namespace std;

namespace OpenMS
{
  const std::string MSQuantifications::NamesOfQuantTypes[] = {"MS1LABEL", "MS2LABEL", "LABELFREE"};

  /// Detailed Constructor
  MSQuantifications::MSQuantifications(const FeatureMap& fm, ExperimentalSettings& es, std::vector<DataProcessing>& dps, std::vector<std::vector<std::pair<String, double> > > label) :
    ExperimentalSettings()
  {
    MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::LABELFREE;
    this->setAnalysisSummaryQuantType(quant_type);

    //~ AssayList,InputFiles,SoftwareList
    //~ aus exp.
    this->registerExperiment(es,dps,std::move(label));
    
    this->setDataProcessingList(fm.getDataProcessing()); //TODO add dp from experiment (i.e. mzml) ?
    feature_maps_  = std::vector<FeatureMap > (1,fm);
  }

  MSQuantifications::~MSQuantifications() = default;

  /// Equality operator
  bool MSQuantifications::operator==(const MSQuantifications & rhs) const
  {
    return ExperimentalSettings::operator==(rhs);
  }

  /// Equality operator
  bool MSQuantifications::operator!=(const MSQuantifications & rhs) const
  {
    return !(operator==(rhs));
  }

  void MSQuantifications::setDataProcessingList(const std::vector<DataProcessing> & dpl)
  {
    data_processings_ = dpl;
  }

  const std::vector<DataProcessing> MSQuantifications::getDataProcessingList() const
  {
    std::vector<DataProcessing> list = data_processings_;

    //This is one way street for dataprocessing - it probably wont get mapped back after write out and reading
    for (std::vector<FeatureMap >::const_iterator fit = feature_maps_.begin(); fit != feature_maps_.end(); ++fit)
    {
      list.insert(list.end(), fit->getDataProcessing().begin(), fit->getDataProcessing().end());
    }

    for (std::vector<ConsensusMap>::const_iterator cit = consensus_maps_.begin(); cit != consensus_maps_.end(); ++cit)
    {
      list.insert(list.end(), cit->getDataProcessing().begin(), cit->getDataProcessing().end());
    }

    return list;
  }

  const std::vector<MSQuantifications::Assay> & MSQuantifications::getAssays() const
  {
    return assays_;
  }

  std::vector<MSQuantifications::Assay> & MSQuantifications::getAssays()
  {
    return assays_;
  }

  //~ std::map<String,ConsensusFeature::Ratio>& MSQuantifications::getRatioCalculations()
  //~ {
  //~ return ratio_calculations_;
  //~ }

  const std::vector<FeatureMap > & MSQuantifications::getFeatureMaps() const
  {
    return feature_maps_;
  }

  const std::vector<ConsensusMap> & MSQuantifications::getConsensusMaps() const
  {
    return consensus_maps_;
  }

  void MSQuantifications::setConsensusMaps(const std::vector<ConsensusMap> & consensus_maps)
  {
    consensus_maps_ = consensus_maps;
  }

  std::vector<ConsensusMap> & MSQuantifications::getConsensusMaps()
  {
    return consensus_maps_;
  }

  const MSQuantifications::AnalysisSummary & MSQuantifications::getAnalysisSummary() const
  {
    return analysis_summary_;
  }

  MSQuantifications::AnalysisSummary & MSQuantifications::getAnalysisSummary()
  {
    return analysis_summary_;
  }

  void MSQuantifications::setAnalysisSummaryQuantType(MSQuantifications::QUANT_TYPES r)
  {
    analysis_summary_.quant_type_ = r;
  }

  void MSQuantifications::addConsensusMap(ConsensusMap & m)
  {
    consensus_maps_.push_back(m);
  }

  void MSQuantifications::assignUIDs()
  {
    for (std::vector<Assay>::iterator ait = assays_.begin(); ait != assays_.end(); ++ait)
    {
      ait->uid_ = String(UniqueIdGenerator::getUniqueId());
    }
  }

  void MSQuantifications::registerExperiment(PeakMap & exp, std::vector<std::vector<std::pair<String, double> > > label)
  {
    for (std::vector<std::vector<std::pair<String, double> > >::const_iterator lit = label.begin(); lit != label.end(); ++lit)
    {
      //TODO look for existing labels
      Assay a;
      a.mods_ = (*lit);
      a.raw_files_.push_back(exp.getExperimentalSettings());
      assays_.push_back(a);
    }

    //TODO check if empty, overwrite MSExperiments inherited front method to work. [0] operator is ugly!
    data_processings_.clear();
    for (Size i = 0; i < exp[0].getDataProcessing().size(); i++)
    {
      data_processings_.push_back( *exp[0].getDataProcessing()[i].get() );
    }
  }
  
  void MSQuantifications::registerExperiment(ExperimentalSettings & es, std::vector<DataProcessing>& /* dps */,  std::vector<std::vector<std::pair<String, double> > > label)
  {
    for (std::vector<std::vector<std::pair<String, double> > >::const_iterator lit = label.begin(); lit != label.end(); ++lit)
    {
      //TODO look for existing labels
      Assay a;
      a.mods_ = (*lit);
      a.raw_files_.push_back(es);
      assays_.push_back(a);
    }
    if (label.empty())
    {
      Assay a;
      a.raw_files_.push_back(es);
      assays_.push_back(a);
    }

    //~ data_processings_ = dps; //TODO add in set fashion
  }

} //namespace OpenMS

