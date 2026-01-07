// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <functional>

namespace OpenMS
{
  /**
      @brief Description of a ion detector (part of a MS Instrument)

      @ingroup Metadata
  */
  class OPENMS_DLLAPI IonDetector :
    public MetaInfoInterface
  {
public:
    /// Detector type
    enum class Type
    {
      TYPENULL,                                  ///< Unknown
      ELECTRONMULTIPLIER,                        ///< Electron multiplier
      PHOTOMULTIPLIER,                           ///< Photo multiplier
      FOCALPLANEARRAY,                           ///< Focal plane array
      FARADAYCUP,                                ///< Faraday cup
      CONVERSIONDYNODEELECTRONMULTIPLIER,        ///< Conversion dynode electron multiplier
      CONVERSIONDYNODEPHOTOMULTIPLIER,           ///< Conversion dynode photo multiplier
      MULTICOLLECTOR,                            ///< Multi-collector
      CHANNELELECTRONMULTIPLIER,                 ///< Channel electron multiplier
      CHANNELTRON,                               ///< channeltron
      DALYDETECTOR,                              ///< daly detector
      MICROCHANNELPLATEDETECTOR,                 ///< microchannel plate detector
      ARRAYDETECTOR,                             ///< array detector
      CONVERSIONDYNODE,                          ///< conversion dynode
      DYNODE,                                    ///< dynode
      FOCALPLANECOLLECTOR,                       ///< focal plane collector
      IONTOPHOTONDETECTOR,                       ///< ion-to-photon detector
      POINTCOLLECTOR,                            ///< point collector
      POSTACCELERATIONDETECTOR,                  ///< postacceleration detector
      PHOTODIODEARRAYDETECTOR,                   ///< photodiode array detector
      INDUCTIVEDETECTOR,                         ///< inductive detector
      ELECTRONMULTIPLIERTUBE,                    ///< electron multiplier tube
      SIZE_OF_TYPE
    };
    /// Names of detector types
    static const std::string NamesOfType[static_cast<size_t>(Type::SIZE_OF_TYPE)];

    /// Acquisition mode
    enum class AcquisitionMode
    {
      ACQMODENULL,             ///< Unknown
      PULSECOUNTING,           ///< Pulse counting
      ADC,                     ///< Analog-digital converter
      TDC,                     ///< Time-digital converter
      TRANSIENTRECORDER,       ///< Transient recorder
      SIZE_OF_ACQUISITIONMODE
    };
    /// Names of acquisition modes
    static const std::string NamesOfAcquisitionMode[static_cast<size_t>(AcquisitionMode::SIZE_OF_ACQUISITIONMODE)];

    /**
     @brief Returns all detector type names known to OpenMS

     @return List of all detector type names
    */
    static StringList getAllNamesOfType();

    /**
     @brief Returns all acquisition mode names known to OpenMS

     @return List of all acquisition mode names
    */
    static StringList getAllNamesOfAcquisitionMode();

    /**
     @brief Convert a Type enum to its string representation

     @param type The detector type enum value to convert
     @return Reference to the string representation
     @throws Exception::InvalidValue if @p type is SIZE_OF_TYPE
    */
    static const std::string& typeToString(Type type);

    /**
     @brief Convert a string to a Type enum

     @param name The string name to convert
     @return The corresponding Type enum value
     @throws Exception::InvalidValue if @p name is not found in NamesOfType[]
    */
    static Type toType(const std::string& name);

    /**
     @brief Convert an AcquisitionMode enum to its string representation

     @param mode The acquisition mode enum value to convert
     @return Reference to the string representation
     @throws Exception::InvalidValue if @p mode is SIZE_OF_ACQUISITIONMODE
    */
    static const std::string& acquisitionModeToString(AcquisitionMode mode);

    /**
     @brief Convert a string to an AcquisitionMode enum

     @param name The string name to convert
     @return The corresponding AcquisitionMode enum value
     @throws Exception::InvalidValue if @p name is not found in NamesOfAcquisitionMode[]
    */
    static AcquisitionMode toAcquisitionMode(const std::string& name);

    /// Constructor
    IonDetector();
    /// Copy constructor
    IonDetector(const IonDetector &) = default;
    /// Move constructor
    IonDetector(IonDetector&&) = default;
    /// Destructor
    ~IonDetector();

    /// Assignment operator
    IonDetector & operator=(const IonDetector &) = default;
    /// Move assignment operator
    IonDetector& operator=(IonDetector&&) & = default;

    /// Equality operator
    bool operator==(const IonDetector & rhs) const;
    /// Equality operator
    bool operator!=(const IonDetector & rhs) const;

    /// returns the detector type
    Type getType() const;
    /// sets the detector type
    void setType(Type type);

    /// returns the acquisition mode
    AcquisitionMode getAcquisitionMode() const;
    /// sets the acquisition mode
    void setAcquisitionMode(AcquisitionMode acquisition_mode);

    /// returns the resolution (in ns)
    double getResolution() const;
    /// sets the resolution (in ns)
    void setResolution(double resolution);

    /// returns the analog-to-digital converter sampling frequency (in Hz)
    double getADCSamplingFrequency() const;
    /// sets the analog-to-digital converter sampling frequency (in Hz)
    void setADCSamplingFrequency(double ADC_sampling_frequency);

    /**
        @brief returns the position of this part in the whole Instrument.

        Order can be ignored, as long the instrument has this default setup:
        - one ion source
        - one or many mass analyzers
        - one ion detector

        For more complex instruments, the order should be defined.
    */
    Int getOrder() const;
    /// sets the order
    void setOrder(Int order);

protected:
    Type type_;
    AcquisitionMode acquisition_mode_;
    double resolution_;
    double ADC_sampling_frequency_;
    Int order_;

  };
} // namespace OpenMS

namespace std
{
  /**
   * @brief Hash function for OpenMS::IonDetector.
   *
   * Hashes type, acquisition mode, resolution, ADC sampling frequency, and order.
   * Note: MetaInfoInterface base class fields are not included in the hash as
   * MetaInfoInterface does not yet have a hash implementation. This hash is
   * consistent with operator== for objects that differ only in the IonDetector-specific
   * fields.
   */
  template<>
  struct hash<OpenMS::IonDetector>
  {
    std::size_t operator()(const OpenMS::IonDetector& d) const noexcept
    {
      std::size_t seed = OpenMS::hash_int(static_cast<int>(d.getType()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(static_cast<int>(d.getAcquisitionMode())));
      OpenMS::hash_combine(seed, OpenMS::hash_float(d.getResolution()));
      OpenMS::hash_combine(seed, OpenMS::hash_float(d.getADCSamplingFrequency()));
      OpenMS::hash_combine(seed, OpenMS::hash_int(d.getOrder()));
      return seed;
    }
  };
} // namespace std

