// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Axel Walter $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <fstream>

namespace OpenMS
{

  class OPENMS_DLLAPI SiriusMSFile
  {
public:

  ///< class to store information about accessions
  class AccessionInfo
  {
  public:
    String sf_path; ///< sourcefile path for mztab-m
    String sf_type; ///< sourcefile type for mztab-m
    String sf_filename; ///< sourcefile name for mztab-m
    String sf_accession; ///< sourcefile accessions for mztab-m
    String native_id_accession; ///< nativeID accession for mztab-m
    String native_id_type; ///< nativeID type for mztab-m
  };

  ///< class to store the compound information
  ///< needed for the mapping of compound and fragment annotated spectrum
  class CompoundInfo
  {
  public:
    String cmp; ///< query_id used compound in .ms file
    double pmass; ///< parent/precursor mass of the compound
    double pint_mono; ///< parent/precursor intensity of the compound
    double rt; ///< retention time of the compound
    double fmz; ///< annotated mass of a feature (if available)
    String fid; ///< annotated feature_id (if available)
    String formula; ///< sumformula of the compound
    int charge; ///< precursor/feature charge
    String ionization; ///< adduct information
    String des; ///< description/name of the compound
    String specref_format; ///< spectra ref format for mztab-m
    String source_file; ///< sourcefile for mztab-m
    String source_format; ///< format of the sourcefile for mztab-m
    std::vector<String> native_ids; ///< native ids of the associated spectra
    String native_ids_id; ///< concatenated list of the associated spectra
    std::vector<String> m_ids; ///< native ids and identifier for multiple possible identification via AMS ("|" separator)
    String m_ids_id; ///< concatenated list of native ids and identifier for multiple possible identification via AMS ("|" separator) used for mapping of compounds and the annotated spectrum.
    std::vector<String> scan_indices; ///< index of the associated spectra
    std::vector<String> specrefs; ///< spectra reference for mztab-m
    int file_index; ///< source file index >
  };

  /**
    @brief Internal structure used in @ref TOPP_SiriusExport that is used
    for the conversion of a MzMlFile to an internal format.

    @ingroup Analysis_ID

    Write content of one mzML/featureXML(optional) file pair to SIRIUS .ms file ofstream.
    Comments (see CompoundInfo) are written to SIRIUS .ms file and additionally stored in CompoundInfo struct.
    If adduct information for a spectrum is missing, no adduct information is added. 
    In this case, SIRIUS assumes default adducts for the respective spectrum.

    @param spectra: Peakmap from input mzML.
    @param os: Write output for .ms file to ofstream.
    @param feature_mapping: Adducts and features (index).
    @param feature_only: Only use features.
    @param isotope_pattern_iterations: At which depth to stop isotope_pattern extraction (if possible).
    @param no_masstrace_info_isotope_pattern: bool if isotope pattern should be extracted (if not in feature)
    @param v_cmpinfo: Vector of CompoundInfo.
    @param file_index: file index (to differentiate entries derived from different mzML files and resolve ambiguities)
    */

    static void store(const MSExperiment& spectra,
                      std::ofstream& os,
                      const FeatureMapping::FeatureToMs2Indices& feature_mapping,
                      const bool& feature_only,
                      const int& isotope_pattern_iterations,
                      const bool no_masstrace_info_isotope_pattern,
                      std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                      const size_t& file_index);
    /**
    @brief Store CompoundInfo objects in tsv file format

    @param v_cmpinfo: Vector with CompoundInfo objects
    @param filename: Filename for tsv file
    */
    static void saveFeatureCompoundInfoAsTSV(const std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                                      const std::string& filename);

  protected:
    /**
    @brief Internal structure to write the .ms file (called in store function)

    @param os: stream
    @param spectra: spectra
    @param ms2_spectra_index: vector of index ms2 spectra (in feature)
    @param ainfo: accession information
    @param adducts: vector of adducts
    @param v_description: vector of descriptions
    @param v_sumformula: vector of sumformulas
    @param f_isotopes: isotope pattern of the feature
    @param feature_charge: feature charge
    @param feature_id: feature id
    @param feature_rt: features retention time
    @param feature_mz: feature mass to charge
    @param writecompound: bool if new compound should be written in .ms file
    @param no_masstrace_info_isotope_pattern: bool if isotope pattern should be extracted (if not in feature)
    @param isotope_pattern_iterations: number of iterations (trying to find a C13 pattern)
    @param count_skipped_spectra: count number of skipped spectra
    @param count_assume_mono: count number of features where mono charge was assumed
    @param count_no_ms1: count number of compounds without a valid ms1 spectrum
    @param v_cmpinfo: vector of CompoundInfo
    @param file_index: file index (to differentiate entries derived from different mzML files and resolve ambiguities)
    */

    static void writeMsFile_(std::ofstream& os,
                             const MSExperiment& spectra,
                             const std::vector<size_t>& ms2_spectra_index,
                             const SiriusMSFile::AccessionInfo& ainfo,
                             const StringList& adducts,
                             const std::vector<String>& v_description,
                             const std::vector<String>& v_sumformula,
                             const std::vector<std::pair<double,double>>& f_isotopes,
                             int& feature_charge,
                             uint64_t& feature_id,
                             const double& feature_rt,
                             const double& feature_mz,
                             bool& writecompound,
                             const bool& no_masstrace_info_isotope_pattern,
                             const int& isotope_pattern_iterations,
                             int& count_skipped_spectra,
                             int& count_assume_mono,
                             int& count_no_ms1,
                             std::vector<SiriusMSFile::CompoundInfo>& v_cmpinfo,
                             const size_t& file_index);
    /**
    @brief Find highest intensity peak near target mz to test if within a margin of error

    @param test_mz: Mass-to-charge to test
    @param spectrum: Spectrum to test
    @param tolerance: Tolerance window (e.g. 10)
    @param ppm: Unit of tolerance window either ppm or Da
    */
    static Int getHighestIntensityPeakInMZRange_(double test_mz,
                                                 const MSSpectrum& spectrum,
                                                 double tolerance,
                                                 bool ppm);

    /**
    @brief Extract precursor isotope pattern if no feature information is available
    based on C12C13 distance.

    @param precursor_mz: Precursor mass-to-charge
    @param precursor_spectrum: Precursor spectrum
    @param iterations: Number of isotopes, which are tried to be extracted
    @param charge: Charge of the precursor
    */
    static std::vector<Peak1D> extractPrecursorIsotopePattern_(const double& precursor_mz,
                                                               const MSSpectrum& precursor_spectrum,
                                                               int& iterations,
                                                               const int& charge);



  };

} // namespace OpenMS