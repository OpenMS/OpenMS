// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Martin Langwisch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <map>
#include <vector>
#include <cmath>

namespace OpenMS
{
  class ProteinIdentification;
  class DateTime;

  /**
  @brief Representation of a Sequest output file

  This class serves to read in a Sequest outfile. The information can be
  retrieved via the load function.

  @todo Handle Modifications (Andreas)
      @todo Complete rewrite of the parser (and those of InsPecT and PepNovo), the code is bullshit... (Andreas)

  @ingroup FileIO
*/
  class OPENMS_DLLAPI SequestOutfile
  {
public:
    /// Constructor
    SequestOutfile();

    /// copy constructor
    SequestOutfile(const SequestOutfile & sequest_outfile);

    /// destructor
    virtual ~SequestOutfile();

    /// assignment operator
    SequestOutfile & operator=(const SequestOutfile & sequest_outfile);

    /// equality operator
    bool operator==(const SequestOutfile & sequest_outfile) const;

    /**
       @brief loads data from a Sequest outfile

       @param result_filename the file to be loaded
       @param peptide_identifications the identifications
       @param protein_identification the protein identifications
       @param p_value_threshold the significance level (for the peptide hit scores)
       @param pvalues a list with the pvalues of the peptides (pvalues computed with peptide prophet)
       @param database the database used for the search
       @param ignore_proteins_per_peptide this is a hack to deal with files that use a suffix like "+1" in column "Reference", but do not actually list extra protein references in subsequent lines

       @throw Exception::FileNotFound is thrown if the given result file could not be found
       @throw Exception::ParseError is thrown if the given result file could not be parsed
       @throw Exception::IllegalArgument

       This class serves to read in a Sequest outfile. The information can be
       retrieved via the load function.
   */
    void load(const String & result_filename, std::vector<PeptideIdentification> & peptide_identifications, ProteinIdentification & protein_identification, const double p_value_threshold, std::vector<double> & pvalues, const String & database = "", const bool ignore_proteins_per_peptide = false);

// /// retrieve the p-values from the out files
//          void getPValuesFromOutFiles(vector< pair < String, vector< double > > >& filenames_and_pvalues) throw (Exception::FileNotFound&, Exception::ParseError);

    /// retrieve columns from a Sequest outfile line
    bool getColumns(const String & line, std::vector<String> & substrings, Size number_of_columns, Size reference_column);

    /** retrieve sequences from a FASTA database
            @param database_filename
            @param ac_position_map
            @param sequences
            @param found
            @param not_found
            @throw Exception::FileNotFound is thrown if the database file could not be found
    */
    void getSequences(const String & database_filename, const std::map<String, Size> & ac_position_map, std::vector<String> & sequences, std::vector<std::pair<String, Size> > & found, std::map<String, Size> & not_found);

    /// retrieve the accession type and accession number from a protein description line
    /// (e.g. from FASTA line: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus], get ac:AAD44166.1 ac type: GenBank)
    void getACAndACType(String line, String & accession, String & accession_type);

    /** read the header of an out file and retrieve various information

            @throw Exception::FileNotFound is thrown if the results file could not be found
            @throw Exception::ParseError is thrown if the results file could not be parsed
    */
    void readOutHeader(const String & result_filename, DateTime & datetime, double & precursor_mz_value, Int & charge, Size & precursor_mass_type, Size & ion_mass_type, Size & displayed_peptides, String & sequest, String & sequest_version, String & database_type, Int & number_column, Int & rank_sp_column, Int & id_column, Int & mh_column, Int & delta_cn_column, Int & xcorr_column, Int & sp_column, Int & sf_column, Int & ions_column, Int & reference_column, Int & peptide_column, Int & score_column, Size & number_of_columns);

private:

    static double const_weights_[];
    static double xcorr_weights_[];
    static double delta_cn_weights_[];
    static double rank_sp_weights_[];
    static double delta_mass_weights_[];
    static Size max_pep_lens_[];
    static Size num_frags_[];
  };

} //namespace OpenMS

