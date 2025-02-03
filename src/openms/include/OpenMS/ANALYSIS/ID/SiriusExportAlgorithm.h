// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka, Lukas Zimmermann $
// --------------------------------------------------------------------------

#pragma once 

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QStringList>

using namespace std;

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusExportAlgorithm : public DefaultParamHandler
    {
    public:
      /// default constructor
      SiriusExportAlgorithm();
      
      // accessor for preprocessing parameters
      bool isFeatureOnly() const { return param_.getValue("feature_only").toBool(); }
      UInt getFilterByNumMassTraces() const { return param_.getValue("filter_by_num_masstraces"); }
      double getPrecursorMzTolerance() const { return param_.getValue("precursor_mz_tolerance"); }
      double getPrecursorRtTolerance() const { return param_.getValue("precursor_rt_tolerance"); }
      bool precursorMzToleranceUnitIsPPM() const { return param_.getValue("precursor_mz_tolerance_unit") == "ppm"; }
      bool isNoMasstraceInfoIsotopePattern() const { return param_.getValue("no_masstrace_info_isotope_pattern").toBool(); }
      int getIsotopePatternIterations() const { return  param_.getValue("isotope_pattern_iterations"); }

      /**
      @brief Preprocessing needed for SIRIUS and AssayGeneratorMetabo

      Filter number of masstraces and perform feature mapping.

      @param featureinfo Path to featureXML
      @param spectra Input of MSExperiment with spectra information
      @param feature_mapping_info Emtpy - stores FeatureMaps and KDTreeMaps internally 
      @param feature_ms2_indices Empty FeatureToMs2Indices
      */
      void preprocessing(const String& featureinfo,
                               const MSExperiment& spectra,
                               FeatureMapping::FeatureMappingInfo& feature_mapping_info,
                               FeatureMapping::FeatureToMs2Indices& feature_ms2_indices) const;

      /**
      @brief logs number of features and spectra used

      Prints the number of features and spectra used (OPENMS_LOG_INFO)

      @param featureinfo Path to featureXML
      @param feature_ms2_indices FeatureToMs2Indices with feature mapping
      @param spectra Input of MSExperiment with spectra information
      */
      void logFeatureSpectraNumber(const String& featureinfo,
                                   const FeatureMapping::FeatureToMs2Indices& feature_ms2_indices,
                                   const MSExperiment& spectra) const;

      /**
      @brief exports SIRIUS .ms file

      Runs SiriusExport with mzML and featureXML (optional) files as input.
      Generates a SIRIUS .ms file and compound info table (optional).

      @param mzML_files List with paths to mzML files
      @param featureXML_files List with paths to featureXML files
      @param out_ms Output file name for SIRIUS .ms file
      @param out_compoundinfo Output file name for tsv file with compound info
      */
      void run(const StringList& mzML_files,
               const StringList& featureXML_files,
               const String& out_ms,
               const String& out_compoundinfo) const;

    };
} // namespace OpenMS
