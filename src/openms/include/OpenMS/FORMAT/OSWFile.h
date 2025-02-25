// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/OSWData.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <array>
#include <map>

namespace OpenMS
{
  /**
    @brief This class serves for reading in and writing OpenSWATH OSW files.

    See OpenSwathOSWWriter for more functionality.

    The reader and writer returns data in a format suitable for PercolatorAdapter.
    OSW files have a flexible data structure. They contain all peptide query
    parameters of TraML/PQP files with the detected and quantified features of
    OpenSwathWorkflow (feature, feature_ms1, feature_ms2 & feature_transition).

    The OSWFile reader extracts the feature information from the OSW file for
    each level (MS1, MS2 & transition) separately and generates Percolator input
    files. For each of the three Percolator reports, OSWFile writer adds a table
    (score_ms1, score_ms2, score_transition) with the respective confidence metrics.
    These tables can be mapped to the corresponding feature tables, are very similar
    to PyProphet results and can thus be used interchangeably.
    
  */
  class OPENMS_DLLAPI OSWFile
  {
  public:
    
    /// query all proteins, not just one with a particular ID
    static constexpr Size ALL_PROTEINS = -1;

    /// opens an OSW file for reading.
    /// @throws Exception::FileNotReadable if @p filename does not exist
    OSWFile(const String& filename);
    OSWFile(const OSWFile& rhs) = default;
    OSWFile& operator=(const OSWFile& rhs) = default;

    /// read data from an SQLLite OSW file into @p swath_result
    /// Depending on the number of proteins, this could take a while. 
    /// @note If you just want the proteins and transitions without peptides and features, use readMinimal().
    void read(OSWData& swath_result);

    /// Reads in transitions and a list of protein names/IDs
    /// but no peptide/feature/transition mapping data (which could be very expensive).
    /// Use in conjunction with on-demand readProtein() to fully populate proteins with peptide/feature data as needed.
    /// @note If you read in all proteins afterwards in one go anyway, using the read() method will be faster (by about 30%)
    void readMinimal(OSWData& swath_result);

    /**
      @brief populates a protein at index @p index  within @p swath_results with Peptides, unless the protein already has peptides

      Internally uses the proteins ID to search for cross referencing peptides and transitions in the OSW file.

      @param swath_result OSWData obtained from the readMinimal() method
      @param index Index into swath_result.getProteins()[index]. Make sure the index is within the vector's size.
      @throws Exception::InvalidValue if the protein at @p index does not have any peptides present in the OSW file
    */
    void readProtein(OSWData& swath_result, const Size index);

    /// for Percolator data read/write operations
    enum class OSWLevel
    {
      MS1,
      MS2,
      TRANSITION,
      SIZE_OF_OSWLEVEL
    };
    static const std::array<std::string, (Size)OSWLevel::SIZE_OF_OSWLEVEL> names_of_oswlevel;

    struct PercolatorFeature
    {
      PercolatorFeature(double score, double qvalue, double pep)
       : score(score), qvalue(qvalue), posterior_error_prob(pep)
       {}
      PercolatorFeature(const PercolatorFeature& rhs) = default;
      
      double score;
      double qvalue;
      double posterior_error_prob;
    };

    /**
      @brief Reads an OSW SQLite file and returns the data on MS1-, MS2- or transition-level
             as ostream (e.g. stringstream or ofstream).
    */
    static void readToPIN(const std::string& filename, const OSWFile::OSWLevel osw_level, std::ostream& pin_output,
                          const double ipf_max_peakgroup_pep, const double ipf_max_transition_isotope_overlap, const double ipf_min_transition_sn);

    /**
    @brief Updates an OpenSWATH OSW SQLite file with the MS1-, MS2- or transition-level results of Percolator.
    */
    static void writeFromPercolator(const std::string& osw_filename, const OSWFile::OSWLevel osw_level, const std::map< std::string, PercolatorFeature >& features);

    /// extract the RUN::ID from the sqMass file
    /// @throws Exception::SqlOperationFailed more than on run exists
    UInt64 getRunID() const;

  protected:
    /** populate transitions of @p swath_result
    
      Clears swath_result entirely (incl. proteins) before adding transitions.
    */
    void readTransitions_(OSWData& swath_result);

    /**
      @brief fill one (@p prot_id) or all proteins into @p swath_result

      @param[out] swath_result Output data. Proteins are cleared before if ALL_PROTEINS is used.
      @param prot_index Using ALL_PROTEINS queries all proteins (could take some time)

    */
    void getFullProteins_(OSWData& swath_result, Size prot_index = ALL_PROTEINS);

    /// set source file and sqMass run-ID
    void readMeta_(OSWData& data);

  private:
    String filename_;       ///< sql file to open/write to
    SqliteConnector conn_;  ///< SQL connection. Stays open as long as this object lives
    bool has_SCOREMS2_;     ///< database contains pyProphet's score_MS2 table with qvalues
  };

} // namespace OpenMS

