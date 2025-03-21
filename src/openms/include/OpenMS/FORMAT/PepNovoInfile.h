// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <map>


namespace OpenMS
{
  /**
      @brief PepNovo input file adapter.

      Creates a PepNovo_PTMs.txt file for PepNovo search.

  @ingroup FileIO
  */
  class OPENMS_DLLAPI PepNovoInfile
  {
public:
    /// default constructor
    PepNovoInfile();

    /// copy constructor
    PepNovoInfile(const PepNovoInfile & pepnovo_infile);

    /// destructor
    virtual ~PepNovoInfile();

    /// assignment operator
    PepNovoInfile & operator=(const PepNovoInfile & pepnovo_infile);

    /// equality operator
    bool operator==(const PepNovoInfile & pepnovo_infile) const;

    /** stores the experiment data in a PepNovo input file that can be used as input for PepNovo shell execution

        @param filename the file which the input file is stored into
        @throw Exception::UnableToCreateFile is thrown if the given file could not be created
    */
    void store(const String & filename);

    /** @brief generates the PepNovo Infile for given fixed and variable modifications			 *
     *
     * @param fixed_mods StringList of fixed modifications unique identifiers
     * @param variable_mods StringList of variable modifications unique identifiers
     */
    void setModifications(const StringList & fixed_mods, const StringList & variable_mods);

    /** @brief return the modifications.
     *
     *  the modification unique identifiers are mapped to the keys used
     *  in the PepNovo Infile (origin+rounded monoisotopic mass of modification ).
     *  (e.g. modification_key_map["K+16"]=="Oxidation (K)" )
     */
    void getModifications(std::map<String, String> & modification_key_map) const;

private:
    ModificationDefinitionsSet mods_;
    std::map<String, String> mods_and_keys_;
    TextFile ptm_file_;


    /** retrieves the name of modification, and generates the corresponding line for the
        PepNovo infile.
        @param modification the modification
        @param variable should be set to true if it variable
   */
    String handlePTMs_(const String & modification, const bool variable);
  };

} // namespace OpenMS

