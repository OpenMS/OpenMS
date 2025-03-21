// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/FORMAT/MSNumpressCoder.h>

#include <vector>

namespace OpenMS
{
  /**
      @brief Options for loading files containing peak data.
  */
  class OPENMS_DLLAPI PeakFileOptions
  {
public:
    ///Default constructor
    PeakFileOptions();
    ///Copy constructor
    PeakFileOptions(const PeakFileOptions &);
    ///Copy assignment
    PeakFileOptions& operator=(const PeakFileOptions&) = default;
    ///Destructor
    ~PeakFileOptions();

    ///@name Meta data and file format option
    //@{
    ///sets whether or not to load only meta data
    void setMetadataOnly(bool only);
    ///returns whether or not to load only meta data
    bool getMetadataOnly() const;

    /// [mzXML only!] Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
    void setForceMQCompatability(bool forceMQ);
    /// [mzXML only!] Whether to write a scan-index and meta data to indicate a Thermo FTMS/ITMS instrument (required to have parameter control in MQ)
    bool getForceMQCompatability() const;

    /// [mzML only!] Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
    void setForceTPPCompatability(bool forceTPP);
    /// [mzML only!] Whether to skip writing the \<isolationWindow\> tag so that TPP finds the correct precursor m/z
    bool getForceTPPCompatability() const;
    //@}

    ///@name Supplemental data option
    //@{
    ///sets whether or not to write supplemental peak data in MzData files
    void setWriteSupplementalData(bool write);
    ///returns whether or not to write supplemental peak data in MzData files
    bool getWriteSupplementalData() const;
    //@}

    ///@name RT range option
    //@{
    ///restricts the range of RT values for peaks to load
    void setRTRange(const DRange<1> & range);
    ///returns @c true if an RT range has been set
    bool hasRTRange() const;
    ///returns the RT range
    const DRange<1> & getRTRange() const;
    //@}

    ///@name m/z range option
    //@{
    ///restricts the range of MZ values for peaks to load
    void setMZRange(const DRange<1> & range);
    ///returns @c true if an MZ range has been set
    bool hasMZRange() const;
    ///returns the MZ range
    const DRange<1> & getMZRange() const;
    //@}

    ///@name Intensity range option
    //@{
    ///restricts the range of intensity values for peaks to load
    void setIntensityRange(const DRange<1> & range);
    ///returns @c true if an intensity range has been set
    bool hasIntensityRange() const;
    ///returns the intensity range
    const DRange<1> & getIntensityRange() const;
    //@}

    /**
        @name MS levels option

        With this option, MS level filters can be set.

        @note The original spectrum identifiers are stored as the nativeID of the spectrum.
    */
    //@{
    ///sets the desired MS levels for peaks to load
    void setMSLevels(const std::vector<Int> & levels);
    ///adds a desired MS level for peaks to load
    void addMSLevel(int level);
    ///clears the MS levels
    void clearMSLevels();
    ///returns @c true, if MS levels have been set
    bool hasMSLevels() const;
    ///returns @c true, if MS level @p level has been set
    bool containsMSLevel(int level) const;
    ///returns the set MS levels
    const std::vector<Int> & getMSLevels() const;
    //@}

    /**
        @name Compression options

        @note This option is ignored if the format does not support compression
    */
    //@{
    //Sets if data should be compressed when writing
    void setCompression(bool compress);
    /// returns @c true, if data should be compressed when writing
    bool getCompression() const;
    //@}

    ///@name lazyload option
    ///sets whether or not to always append the data to the given map (even if a consumer is given)
    void setAlwaysAppendData(bool only);
    ///returns whether or not to always append the data to the given map (even if a consumer is given)
    bool getAlwaysAppendData() const;
    ///sets whether to fill the actual data into the container (spectrum/chromatogram)
    void setFillData(bool only);
    ///returns whether to fill the actual data into the container (spectrum/chromatogram)
    bool getFillData() const;
    ///sets whether to skip some XML checks and be fast instead
    void setSkipXMLChecks(bool only);
    ///returns whether to skip some XML checks and be fast instead
    bool getSkipXMLChecks() const;

    /// @name sort peaks in spectra / chromatograms by position
    ///sets whether or not to sort peaks in spectra
    void setSortSpectraByMZ(bool sort);
    ///gets whether or not peaks in spectra should be sorted
    bool getSortSpectraByMZ() const;
    ///sets whether or not to sort peaks in chromatograms
    void setSortChromatogramsByRT(bool sort);
    ///gets whether or not peaks in chromatograms should be sorted
    bool getSortChromatogramsByRT() const;

    /**
        @name Precision options

        Note that m/z and RT are controlled with the same flag (for spectra and
        chromatograms) while there is a separate flag for intensity.

        @note This option is ignored if the format does not support multiple precisions
    */
    //@{
    //Sets if mz-data and rt-data should be stored with 32bit or 64bit precision
    void setMz32Bit(bool mz_32_bit);
    //returns @c true, if mz-data and rt-data should be stored with 32bit precision
    bool getMz32Bit() const;
    //Sets if intensity data should be stored with 32bit or 64bit precision
    void setIntensity32Bit(bool int_32_bit);
    //returns @c true, if intensity data should be stored with 32bit precision
    bool getIntensity32Bit() const;
    //@}

    /// Whether to write an index at the end of the file (e.g. indexedmzML file format)
    bool getWriteIndex() const;
    /// Whether to write an index at the end of the file (e.g. indexedmzML file format)
    void setWriteIndex(bool write_index);

    /// Set numpress configuration options for m/z or rt dimension
    MSNumpressCoder::NumpressConfig getNumpressConfigurationMassTime() const;
    /// Get numpress configuration options for m/z or rt dimension
    void setNumpressConfigurationMassTime(MSNumpressCoder::NumpressConfig config);
    /// Set numpress configuration options for intensity dimension
    MSNumpressCoder::NumpressConfig getNumpressConfigurationIntensity() const;
    /// Get numpress configuration options for intensity dimension
    void setNumpressConfigurationIntensity(MSNumpressCoder::NumpressConfig config);
    /// Set numpress configuration options for float data arrays
    MSNumpressCoder::NumpressConfig getNumpressConfigurationFloatDataArray() const;
    /// Get numpress configuration options for float data arrays
    void setNumpressConfigurationFloatDataArray(MSNumpressCoder::NumpressConfig config);

    /**
        @name Data pool size options

        Some file readers and writers can process the data in parallel by
        reading in parts of the file and keeping it in memory and then process
        this partial data in parallel. This parameter specifies how many
        data points (spectra/chromatograms) should be read before parallel
        processing is initiated.
    */
    //@{
    /// Get maximal size of the data pool
    Size getMaxDataPoolSize() const;
    /// Set maximal size of the data pool
    void setMaxDataPoolSize(Size size);
    //@}

    /// [mzML only!] Whether to use the "selected ion m/z" value as the precursor m/z value (alternative: use the "isolation window target m/z" value)
    bool getPrecursorMZSelectedIon() const;

    /// [mzML only!] Set whether to use the "selected ion m/z" value as the precursor m/z value (alternative: use the "isolation window target m/z" value)
    void setPrecursorMZSelectedIon(bool choice);

    /// do these options skip spectra or chromatograms due to RT or MSLevel filters?
    bool hasFilters() const;

private:
    bool metadata_only_;
    bool force_maxquant_compatibility_; ///< for mzXML-writing only: set a fixed vendor (Thermo Scientific), mass analyzer (FTMS)
    bool force_tpp_compatibility_; ///< for mzML-writing only: work around some bugs in TPP file parsers
    bool write_supplemental_data_;
    bool has_rt_range_;
    bool has_mz_range_;
    bool has_intensity_range_;
    bool mz_32_bit_;
    bool int_32_bit_;
    DRange<1> rt_range_;
    DRange<1> mz_range_;
    DRange<1> intensity_range_;
    std::vector<Int> ms_levels_;
    bool zlib_compression_;
    bool always_append_data_;
    bool skip_xml_checks_;
    bool sort_spectra_by_mz_;
    bool sort_chromatograms_by_rt_;
    bool fill_data_;
    bool write_index_;
    MSNumpressCoder::NumpressConfig np_config_mz_;
    MSNumpressCoder::NumpressConfig np_config_int_;
    MSNumpressCoder::NumpressConfig np_config_fda_;
    Size maximal_data_pool_size_;
    bool precursor_mz_selected_ion_;
  };

} // namespace OpenMS

