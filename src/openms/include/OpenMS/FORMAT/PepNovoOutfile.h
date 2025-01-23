// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Sandro Andreotti, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
  class ProteinIdentification;

  /**
      @brief Representation of a PepNovo output file

      This class serves to read in a PepNovo outfile. The information can be
      retrieved via the load function.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI PepNovoOutfile
  {
public:
    
    typedef std::map<Size, std::pair<double, double> > IndexPosMappingType;
    
    /// Constructor
    PepNovoOutfile();

    /// copy constructor
    PepNovoOutfile(const PepNovoOutfile & pepnovo_outfile);

    /// destructor
    virtual ~PepNovoOutfile();

    /// assignment operator
    PepNovoOutfile & operator=(const PepNovoOutfile & pepnovo_outfile);

    /// equality operator
    bool operator==(const PepNovoOutfile & pepnovo_outfile) const;

    /**
       @brief loads data from a PepNovo outfile

       @param result_filename the file to be loaded
       @param peptide_identifications the peptide identifications
       @param protein_identification the protein identification
       @param score_threshold cutoff threshold for the PepNovo score (PnvScr)
       @param id_rt_mz map the spectrum identifiers returned by PepNovo
       to the rt and mz values of the spectrum (used to map the identifications back to the spectra). key= &lt;PepNovo Id&gt;, value= &lt;pair&lt;rt,mz&gt; &gt;.
       For spectra not present in this map identifications cannot be mapped back.
       @param mod_id_map map the OpenMS id for modifications (FullId) to the ids returned by PepNovo key= &lt;PepNovo_key&gt;, value= &lt;OpenMS FullId&gt;
   */
    void load(const std::string & result_filename, std::vector<PeptideIdentification> & peptide_identifications,
              ProteinIdentification & protein_identification,
              const double & score_threshold,
              const IndexPosMappingType & id_rt_mz,
              const std::map<String, String> & mod_id_map);

    /** @brief get the search engine version and search parameters from a PepNovo output file
     *
     * search parameters (precursor tolerance, peak mass tolerance, allowed modifications)are stored in the protein_identification.

        @param pepnovo_output_without_parameters_filename
        @param protein_identification
    */
    void getSearchEngineAndVersion(const String & pepnovo_output_without_parameters_filename, ProteinIdentification & protein_identification);

  };

} //namespace OpenMS

