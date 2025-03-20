// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS
{
  /**
    @brief XTandem input file.

    This class is able to load/write a X!Tandem configuration file.

    These files store parameters within 'note' tags, e.g.,
    @verbatim
      <note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
    	<note type="input" label="output, proteins">yes</note>
    @endverbatim

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XTandemInfile :
    public Internal::XMLFile
  {
public:

    /// error unit, either Da or ppm
    enum ErrorUnit
    {
      DALTONS = 0,
      PPM
    };

    /// Mass type of the precursor, either monoisotopic or average
    enum MassType
    {
      MONOISOTOPIC = 0,
      AVERAGE
    };

    /// constructor
    XTandemInfile();

    /// constructor
    ~XTandemInfile() override;

    /// setter for the fragment mass tolerance
    void setFragmentMassTolerance(double tolerance);

    /// returns the fragment mass tolerance
    double getFragmentMassTolerance() const;

    /// sets the precursor mass tolerance (plus only)
    void setPrecursorMassTolerancePlus(double tol);

    /// returns the precursor mass tolerance (plus only)
    double getPrecursorMassTolerancePlus() const;

    /// set the precursor mass tolerance (minus only)
    void setPrecursorMassToleranceMinus(double tol);

    /// returns the precursor mass tolerance (minus only)
    double getPrecursorMassToleranceMinus() const;

    /// sets the precursor mass type
    void setPrecursorErrorType(MassType mono_isotopic);

    /// returns the precursor mass type
    MassType getPrecursorErrorType() const;

    /// sets the fragment mass error unit (Da, ppm)
    void setFragmentMassErrorUnit(ErrorUnit unit);

    /// returns the fragment mass error unit (Da, ppm)
    ErrorUnit getFragmentMassErrorUnit() const;

    /// sets the precursor mass error unit (Da, ppm)
    void setPrecursorMassErrorUnit(ErrorUnit unit);

    /// returns the precursor mass error unit (Da, ppm)
    ErrorUnit getPrecursorMassErrorUnit() const;

    /// sets the number of threads used during the identifications
    void setNumberOfThreads(UInt threads);

    /// returns the number of threads
    UInt getNumberOfThreads() const;

    /// sets the modifications using a modification definitions set
    void setModifications(const ModificationDefinitionsSet& mods);

    /// returns the modifications set, using a modification definitions set
    const ModificationDefinitionsSet& getModifications() const;

    /// sets the output filename
    void setOutputFilename(const String& output);

    /// returns the output filename
    const String& getOutputFilename() const;

    /// sets the input filename
    void setInputFilename(const String& input_file);

    /// returns the input filename
    const String& getInputFilename() const;

    /// set the filename of the taxonomy file
    void setTaxonomyFilename(const String& filename);

    /// returns the filename of the taxonomy file
    const String& getTaxonomyFilename() const;

    /// sets the default parameters file
    void setDefaultParametersFilename(const String& filename);

    /// returns the default parameters file
    const String& getDefaultParametersFilename() const;

    /// sets the taxon used in the taxonomy file
    void setTaxon(const String& taxon);

    /// returns the taxon used in the taxonomy file
    const String& getTaxon() const;

    /// sets the max precursor charge
    void setMaxPrecursorCharge(Int max_charge);

    /// returns the max precursor charge
    Int getMaxPrecursorCharge() const;

    /// sets the number of missed cleavages allowed
    void setNumberOfMissedCleavages(UInt missed_cleavages);

    /// returns the number of missed cleavages allowed
    UInt getNumberOfMissedCleavages() const;

    /// sets the output result type ("all", "valid" or "stochastic")
    void setOutputResults(const String& result);

    /// returns the output result type ("all", "valid" or "stochastic")
    String getOutputResults() const;

    /// sets the max valid E-value allowed in the list
    void setMaxValidEValue(double value);

    /// returns the max valid E-value allowed in the list
    double getMaxValidEValue() const;

    /// set state of semi cleavage
    void setSemiCleavage(const bool semi_cleavage);

    /// set if misassignment of precursor to first and second 13C isotopic peak should also be considered
    void setAllowIsotopeError(const bool allow_isotope_error);

    /// get state of noise suppression
    bool getNoiseSuppression() const;

    /// set state of noise suppression
    void setNoiseSuppression(const bool noise_suppression);

    /// set the cleavage site with a X! Tandem conform regex
    void setCleavageSite(const String& cleavage_site);
    
    /// returns the cleavage site regex
    const String& getCleavageSite() const;

    /** 
      @brief Writes the X! Tandem input file to the given filename

      If @p ignore_member_parameters is true, only a very limited number of
      tags fed by member variables (i.e. in, out, database/taxonomy) is written.
      
      @param filename the name of the file which is written
      @param ignore_member_parameters Do not write tags for class members
      @param force_default_mods Force writing of mods covered by special parameters
      @throw UnableToCreateFile is thrown if the given file could not be created
    */
    void write(const String& filename, bool ignore_member_parameters = false,
               bool force_default_mods = false);

protected:

    XTandemInfile(const XTandemInfile& rhs);

    XTandemInfile& operator=(const XTandemInfile& rhs);

    void writeTo_(std::ostream& os, bool ignore_member_parameters);

    void writeNote_(std::ostream& os, const String& label, const String& value);

    void writeNote_(std::ostream& os, const String& label, const char* value);

    void writeNote_(std::ostream& os, const String& label, bool value);

    /**
      @brief Converts the given set of Modifications into a format compatible to X!Tandem.

      The set affected_origins can be used to avoid duplicate modifications, which are not supported in X! Tandem.
      Currently, a warning message is printed.
      Also, if a fixed mod is already given, a corresponding variable mods needs to have its delta mass reduced by the fixed modifications mass.
      This is also done automatically here.

      @param mods The modifications to convert
      @param affected_origins Set of origins, which were used previously. Will be augmented with the current mods.

      @return An X! Tandem compatible string representation.
    */
    String convertModificationSet_(const std::set<ModificationDefinition>& mods, std::map<String, double>& affected_origins) const;

    double fragment_mass_tolerance_;

    double precursor_mass_tolerance_plus_;

    double precursor_mass_tolerance_minus_;

    ErrorUnit fragment_mass_error_unit_;

    ErrorUnit precursor_mass_error_unit_;

    MassType fragment_mass_type_;

    MassType precursor_mass_type_;

    UInt max_precursor_charge_;

    double precursor_lower_mz_;

    double fragment_lower_mz_;

    UInt number_of_threads_;

    UInt batch_size_;

    ModificationDefinitionsSet modifications_;

    String input_filename_;

    String output_filename_;

    String taxonomy_file_;

    String taxon_;

    String cleavage_site_;

    /// semi cleavage
    bool semi_cleavage_;

    bool allow_isotope_error_;

    // scoring
    UInt number_of_missed_cleavages_;

    String default_parameters_file_;

    // output parameters
    String output_results_;

    double max_valid_evalue_;

    // force writing of mods covered by special parameters?
    bool force_default_mods_;

  };

} // namespace OpenMS

