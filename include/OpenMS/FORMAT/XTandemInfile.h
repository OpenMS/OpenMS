// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_XTANDEMINFILE_H
#define OPENMS_FORMAT_XTANDEMINFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/HANDLERS/XTandemInfileXMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/FORMAT/XMLFile.h>

namespace OpenMS
{
  /**
      @brief XTandem input file adapter

      This class is able to create a X!Tandem configuration file for a search

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
    virtual ~XTandemInfile();

    //<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
    //<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
    //<note type="input" label="spectrum, parent monoisotopic mass error plus">100</note>
    //<note type="input" label="spectrum, parent monoisotopic mass error minus">100</note>
    //<note type="input" label="spectrum, parent monoisotopic mass isotope error">yes</note>
    //<note type="input" label="spectrum, fragment monoisotopic mass error units">Daltons</note>
    //<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
    //<note type="input" label="spectrum, parent monoisotopic mass error units">ppm</note>
    //<note>The value for this parameter may be 'Daltons' or 'ppm': all other values are ignored</note>
    //<note type="input" label="spectrum, fragment mass type">monoisotopic</note>
    //<note>values are monoisotopic|average </note>

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
    void setModifications(const ModificationDefinitionsSet & mods);

    /// returns the modifications set, using a modification definitions set
    const ModificationDefinitionsSet & getModifications() const;

    /// sets the output filename
    void setOutputFilename(const String & output);

    /// returns the output filename
    const String & getOutputFilename() const;

    /// sets the input filename
    void setInputFilename(const String & input_file);

    /// returns the input filename
    const String & getInputFilename() const;

    /// set the filename of the taxonomy file
    void setTaxonomyFilename(const String & filename);

    /// returns the filename of the taxonomy file
    const String & getTaxonomyFilename() const;

    /// sets the default paramters file
    void setDefaultParametersFilename(const String & filename);

    /// returns the default parameters file
    const String & getDefaultParametersFilename() const;

    /// sets the taxon used in the taxonomy file
    void setTaxon(const String & taxon);

    /// returns the taxon used in the taxonomy file
    const String & getTaxon() const;

    /// sets the max precursor charge
    void setMaxPrecursorCharge(Int max_charge);

    /// returns the max precursor charge
    Int getMaxPrecursorCharge() const;

    /// sets the number of missed cleavages allowed
    void setNumberOfMissedCleavages(UInt missed_cleavages);

    /// returns the number of missed cleavages allowed
    UInt getNumberOfMissedCleavages() const;

    /// sets the max valid E-value allowed in the list
    void setMaxValidEValue(double value);

    /// returns the max valid E-value allowed in the list
    double getMaxValidEValue() const;

    /// get state of refine setting
    bool isRefining() const;

    /// set state of semi cleavage
    void setSemiCleavage(const bool semi_cleavage);

    /// set state of refine setting
    void setRefine(const bool refine);

    /** writes the XTandemInfile to the given file

            @param filename the name of the file which is written
            @throw UnableToCreateFile is thrown if the given file could not be created
    */
    void write(const String & filename);

    /** read the information from the given filename

            @param filename the file which should be read from
            @throw FileNotFound is thrown if the given file could not be found
            @throw ParseError is thrown if the given file could not be parsed
    */
    void load(const String & filename);

protected:

    XTandemInfile(const XTandemInfile & rhs);

    XTandemInfile & operator=(const XTandemInfile & rhs);

    void writeTo_(std::ostream & os);

    void writeNote_(std::ostream & os, const String & type, const String & label, const String & value);

    void writeNote_(std::ostream & os, const String & type, const String & label, const char * value);

    void writeNote_(std::ostream & os, const String & type, const String & label, bool value);

    double fragment_mass_tolerance_;

    double precursor_mass_tolerance_plus_;

    double precursor_mass_tolerance_minus_;

    MassType precursor_mass_type_;

    ErrorUnit precursor_mass_error_unit_;

    ErrorUnit fragment_mass_error_unit_;

    MassType fragment_mass_type_;

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

    // Sment
    bool refine_;

    //semi cleavage
    bool semi_cleavage_;

    double refine_max_valid_evalue_;

    // scoring
    UInt number_of_missed_cleavages_;

    String default_parameters_file_;

    // output parameters
    double max_valid_evalue_;

    //<note type="input" label="spectrum, fragment monoisotopic mass error">0.4</note>
    std::vector<Internal::XTandemInfileNote> notes_;
  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_XTANDEMINFILE_H
