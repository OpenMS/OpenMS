// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>


namespace OpenMS
{
  
  /**
    @brief Base class for Search Engine Adapters

    It is build on top of TOPPBase and provides convenience functions for regular tasks in SearchEngines.

    This base class enforces a common parameter scheme upon each adapter.
    E.g. '-database' and '-in'.
    This might be extended/changed in the future.
    
  */
  class OPENMS_DLLAPI SearchEngineBase : public TOPPBase
  {
  public:
    /// No default constructor
    SearchEngineBase() = delete;

    /// No default copy constructor.
    SearchEngineBase(const SearchEngineBase&) = delete;

    /**
      @brief Constructor

      Must match TOPPBase' Ctor!

      @param name Tool name.
      @param description Short description of the tool (one line).
      @param official If this is an official TOPP tool contained in the OpenMS/TOPP release.
             If @em true the tool name is checked against the list of TOPP tools and a warning printed if missing.
      @param citations Add one or more citations if they are associated specifically to this TOPP tool; they will be printed during `--help`
      @param toolhandler_test Check if this tool is registered with the ToolHandler (disable for unit tests only)
    */
    SearchEngineBase(const String& name, const String& description, bool official = true, const std::vector<Citation>& citations = {}, bool toolhandler_test = true);

    /// Destructor
    ~SearchEngineBase() override;

    /**
      @brief Reads the '-in' argument from internal parameters (usually an mzML file) and checks if MS2 spectra are present and are centroided.

      If the file is an mzML file, the spectra annotation can be checked. If no MS2 or profile MS2 data is found, an exception is thrown.
      If the file is any other format, the overhead of reading in the file is too large and we just issue a general warning that centroided data should be used.

      @param ms_level The MS level to check for their type (centroided/profile)

      @return A filename (might be a relative or absolute path)

      @throws OpenMS::Exception::FileEmpty if no spectra are found (mzML only)
      @throws OpenMS::Exception::IllegalArgument if spectra are not centroided (mzML only)

    */
    String getRawfileName(int ms_level = 2) const;


    /**
      @brief Reads the '-database' argument from internal parameters (or from @p db) and tries to find the db in search directories (if it cannot be found immediately). If not found, an exception is thrown.
      
      @param db [Optional] Instead of reading the '-database', you can provide a custom name here (might be required for special db formats, see OMSSA)
      @return filename for DB (might be a relative or absolute path)
 
      @throws OpenMS::Exception::FileNotFound if database name could not be resolved
    */
    String getDBFilename(const String& db = "") const;


    /**
      @brief Adds option to reassociate peptides with proteins (and annotate target/decoy information)

      @param peptide_indexing_parameter peptide indexer settings. May be modified to enable search engine specific defaults (e.g., not-tryptic etc.). 
    */
    virtual void registerPeptideIndexingParameter_(Param peptide_indexing_parameter);

    /**
      @brief Reindex peptide to protein association
    */
    virtual SearchEngineBase::ExitCodes reindex_(
      std::vector<ProteinIdentification>& protein_identifications, 
      std::vector<PeptideIdentification>& peptide_identifications) const;
  }; // end SearchEngineBase

}   // end NS OpenMS
