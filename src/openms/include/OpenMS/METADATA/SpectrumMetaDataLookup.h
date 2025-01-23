// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/SpectrumLookup.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <limits> // for "quiet_NaN"

namespace OpenMS
{

  /**
    @brief Helper class for looking up spectrum meta data

    The class deals with meta data of spectra and provides functions for the extraction and look-up of this data.

    A common use case for this functionality is importing peptide/protein
    identification results from search engine-specific file formats, where some
    meta information may have to be looked up in the raw data (primarily
    retention times).
    One example of this is in the function addMissingRTsToPeptideIDs().

    Meta data of a spectra is stored in SpectrumMetaDataLookup::SpectrumMetaData structures.
    In order to control which meta data to extract/look-up, flags
    (SpectrumMetaDataLookup::MetaDataFlags) are used.
    Meta data can be extracted from spectra or from spectrum reference strings.
    The format of a spectrum reference is defined via a regular expression
    containing named groups (format "(?<GROUP>...)" for the different data
    items.
    The table below illustrates the different meta data types and how they are represented.

    <CENTER>
        <table>
            <tr>
                <td ALIGN="center" BGCOLOR="#EBEBEB"> @p SpectrumMetaData member </td>
                <td ALIGN="center" BGCOLOR="#EBEBEB"> @p MetaDataFlags flag </td>
                <td ALIGN="center" BGCOLOR="#EBEBEB"> Reg. exp. group </td>
                <td ALIGN="center" BGCOLOR="#EBEBEB"> Comment (*: undefined for MS1 spectra)</td>
            </tr>
            <tr>
                <td ALIGN="center"> @p rt </td>
                <td ALIGN="center"> @p MDF_RT </td>
                <td ALIGN="center"> @p RT </td>
                <td ALIGN="center"> Retention time of the spectrum </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p precursor_rt </td>
                <td ALIGN="center"> @p MDF_PRECURSORRT </td>
                <td ALIGN="center"> @p PRECRT </td>
                <td ALIGN="center"> Retention time of the precursor spectrum* </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p precursor_mz </td>
                <td ALIGN="center"> @p MDF_PRECURSORMZ </td>
                <td ALIGN="center"> @p MZ </td>
                <td ALIGN="center"> Mass-to-charge ratio of the precursor ion* </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p precursor_charge </td>
                <td ALIGN="center"> @p MDF_PRECURSORCHARGE </td>
                <td ALIGN="center"> @p CHARGE </td>
                <td ALIGN="center"> Charge of the precursor ion* </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p ms_level </td>
                <td ALIGN="center"> @p MDF_MSLEVEL </td>
                <td ALIGN="center"> @p LEVEL </td>
                <td ALIGN="center"> MS level (1 for survey scan, 2 for fragment scan, etc.) </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p scan_number </td>
                <td ALIGN="center"> @p MDF_SCANNUMBER </td>
                <td ALIGN="center"> @p SCAN </td>
                <td ALIGN="center"> Scan number (extracted from the native ID) </td>
            </tr>
            <tr>
                <td ALIGN="center"> @p native_id </td>
                <td ALIGN="center"> @p MDF_NATIVEID </td>
                <td ALIGN="center"> @p ID </td>
                <td ALIGN="center"> Native ID of the spectrum </td>
            </tr>
            <tr>
                <td ALIGN="center"></td>
                <td ALIGN="center"> @p MDF_ALL </td>
                <td ALIGN="center"></td>
                <td ALIGN="center"> Shortcut for "all flags set" </td>
            </tr>
            <tr>
                <td ALIGN="center"></td>
                <td ALIGN="center"></td>
                <td ALIGN="center"> @p INDEX0 </td>
                <td ALIGN="center"> Only for look-up: index (vector pos.) counting from 0 </td>
            </tr>
            <tr>
                <td ALIGN="center"></td>
                <td ALIGN="center"></td>
                <td ALIGN="center"> @p INDEX1 </td>
                <td ALIGN="center"> Only for look-up: index (vector pos.) counting from 1 </td>
            </tr>
        </table>
    </CENTER>

    @see OpenMS::SpectrumLookup
  */
  class OPENMS_DLLAPI SpectrumMetaDataLookup: public SpectrumLookup
  {
  public:

    /// Bit mask for which meta data to extract from a spectrum
    typedef unsigned char MetaDataFlags;

    /// Possible meta data to extract from a spectrum. 
    /// Note that the static variables need to be put on separate lines due to a compiler bug in VS
    static const MetaDataFlags MDF_RT = 1;
    static const MetaDataFlags MDF_PRECURSORRT = 2;
    static const MetaDataFlags MDF_PRECURSORMZ = 4;
    static const MetaDataFlags MDF_PRECURSORCHARGE = 8;
    static const MetaDataFlags MDF_MSLEVEL = 16;
    static const MetaDataFlags MDF_SCANNUMBER = 32;
    static const MetaDataFlags MDF_NATIVEID = 64;
    static const MetaDataFlags MDF_ALL = 127;


    /// Meta data of a spectrum
    struct SpectrumMetaData
    {
      double rt; ///< Retention time
      double precursor_rt; ///< Precursor retention time
      double precursor_mz; ///< Precursor mass-to-charge ratio
      Int precursor_charge; ///< Precursor charge
      Size ms_level; ///< MS level
      Int scan_number; ///< Scan number
      String native_id; ///< Native ID

      /// Constructor
      SpectrumMetaData():
        rt(std::numeric_limits<double>::quiet_NaN()), 
        precursor_rt(std::numeric_limits<double>::quiet_NaN()), 
        precursor_mz(std::numeric_limits<double>::quiet_NaN()),
        precursor_charge(0), ms_level(0), scan_number(-1), native_id("")
      {
      }

      /// Copy constructor
      SpectrumMetaData(const SpectrumMetaData &) = default;
      /// Move constructor
      SpectrumMetaData(SpectrumMetaData&&) = default;
      /// Destructor
      ~SpectrumMetaData() = default;

      /// Assignment operator
      SpectrumMetaData & operator=(const SpectrumMetaData &) = default;
      /// Move assignment operator
      SpectrumMetaData& operator=(SpectrumMetaData&&) & = default;
    };

    /// Constructor
    SpectrumMetaDataLookup(): SpectrumLookup()
    {}

    /// Destructor
    ~SpectrumMetaDataLookup() override {}

    /**
       @brief Read spectra and store their meta data

       @tparam SpectrumContainer Spectrum container class, must support @p size and @p operator[]
       @param spectra Container of spectra
       @param scan_regexp Regular expression for matching scan numbers in spectrum native IDs (must contain the named group "?<SCAN>")
       @param get_precursor_rt Assign precursor retention times? (This relies on all precursor spectra being present and in the right order.)

       @throw Exception::IllegalArgument if @p scan_regexp does not contain "?<SCAN>" (and is not empty)
    */
    template <typename SpectrumContainer>
    void readSpectra(const SpectrumContainer& spectra, 
                     const String& scan_regexp = default_scan_regexp,
                     bool get_precursor_rt = false)
    {
      // If class SpectrumContainer is e.g. OnDiscMSExperiment, reading each
      // spectrum is expensive. Thus, to avoid iterating over all spectra twice,
      // we do not call "SpectrumLookup::readSpectra" here:
      n_spectra_ = spectra.size();
      metadata_.reserve(n_spectra_);
      setScanRegExp_(scan_regexp);
      // mapping: MS level -> RT of previous spectrum of that level
      std::map<Size, double> precursor_rts;
      for (Size i = 0; i < n_spectra_; ++i)
      {
        const MSSpectrum& spectrum = spectra[i];
        SpectrumMetaData meta;
        getSpectrumMetaData(spectrum, meta, scan_regexp_, precursor_rts);
        if (get_precursor_rt) precursor_rts[meta.ms_level] = meta.rt;
        addEntry_(i, meta.rt, meta.scan_number, meta.native_id);
        metadata_.push_back(meta);
      }
    }


    /**
     * @brief set spectra_data from read SpectrumContainer origin (i.e. filename)
     *
     * @param spectra_data the name (and path) of the origin of the read SpectrumContainer
     */
    void setSpectraDataRef(const String& spectra_data)
    {
      this->spectra_data_ref = spectra_data;
    }


    /**
       @brief Look up meta data of a spectrum

       @param index Index of the spectrum
       @param meta Meta data output
    */
    void getSpectrumMetaData(Size index, SpectrumMetaData& meta) const;

    /**
       @brief Extract meta data from a spectrum

       @param spectrum Spectrum input
       @param meta Meta data output
       @param scan_regexp Regular expression for extracting scan number from spectrum native ID
       @param precursor_rts RTs of potential precursor spectra of different MS levels

       Scan number and precursor RT, respectively, are only extracted if @p scan_regexp/@p precursor_rts are not empty.
    */
    static void getSpectrumMetaData(
      const MSSpectrum& spectrum, SpectrumMetaData& meta,
      const boost::regex& scan_regexp = boost::regex(),
      const std::map<Size, double>& precursor_rts = (std::map<Size, double>()));

    /**
       @brief Extract meta data via a spectrum reference

       @param spectrum_ref Spectrum reference to parse
       @param meta Meta data output
       @param flags What meta data to extract

       @throw Exception::ElementNotFound if a spectrum look-up was necessary, but no matching spectrum was found

       This function is a combination of getSpectrumMetaData() and SpectrumLookup::findByReference(). However, the spectrum is only looked up if necessary, i.e. if the required meta data - as defined by @p flags - cannot be extracted from the spectrum reference itself.
    */
    void getSpectrumMetaData(const String& spectrum_ref, SpectrumMetaData& meta,
                             MetaDataFlags flags = MDF_ALL) const;

	/**
	   @brief Add missing retention time (RT) values to peptide identifications based on raw data

	   @param peptides Peptide IDs with or without RT values
	   @param exp The MSExperiment object representing the raw data file (e.g., mzML) used to look up RT values.

	   @return True if all peptide IDs could be annotated successfully (including if all already had RT values), false otherwise.

	*/
	static bool addMissingRTsToPeptideIDs(std::vector<PeptideIdentification>& peptides, const MSExperiment& exp);

    /**
     * @brief Adds missing ion mobility information to peptide identifications.
     * 
     * This function adds missing ion mobility (IM) information to the peptide identifications.
     * The missing IM information is retrieved from the MSExperiment.
     * 
     * @param peptides The vector of peptide identifications to update.
     * @param exp The MSExperiment object representing the raw data file (e.g., mzML) used to look up IM values.
     * 
     * @return True if all missing IM information was successfully added to the peptide identifications, false otherwise.
    */
    static bool addMissingIMToPeptideIDs(std::vector<PeptideIdentification>& peptides,
    									const MSExperiment& exp);

    /**
     * @brief Add missing "spectrum_reference"s to peptide identifications based on raw data
     *
     * @param peptides Peptide IDs with or without spectrum_reference
     * @param filename the name of the mz_file from which to draw spectrum_references
     * @param stop_on_error Stop when an ID could not be matched to a spectrum (or keep going)?
     * @param override_spectra_data if given ProteinIdentifications should be updated with new "spectra_data" values from SpectrumMetaDataLookup
     * @param override_spectra_references if given PeptideIdentifications with existing spectrum_reference should be updated from SpectrumMetaDataLookup
     * @param proteins Protein IDs corresponding to the Peptide IDs
     *
     * @return True if all peptide IDs could be annotated successfully (including if all already had "spectrum_reference" values), false otherwise.
     *
     * Look-up works by matching RT of a peptide identification with the given spectra. Matched spectra 'native ID' will be annotated to the identification. All spectrum_references are updated/added.
     */
    static bool addMissingSpectrumReferences(std::vector<PeptideIdentification>& peptides, 
      const String& filename,
      bool stop_on_error = false, 
      bool override_spectra_data = false, 
      bool override_spectra_references = false, 
      std::vector<ProteinIdentification> proteins = std::vector<ProteinIdentification>());



  protected:

    std::vector<SpectrumMetaData> metadata_; ///< Meta data for spectra
    String spectra_data_ref;

  private:

    /// Copy constructor (not implemented)
    SpectrumMetaDataLookup(const SpectrumMetaDataLookup&);

    /// Assignment operator (not implemented)
    SpectrumMetaDataLookup& operator=(const SpectrumMetaDataLookup&);

  };

} //namespace OpenMS

