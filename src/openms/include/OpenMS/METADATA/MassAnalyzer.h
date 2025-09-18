// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Description of a mass analyzer (part of a MS Instrument)

      @ingroup Metadata
  */
  class OPENMS_DLLAPI MassAnalyzer :
    public MetaInfoInterface
  {
public:
    /// analyzer type
    enum AnalyzerType
    {
      ANALYZERNULL,                 ///< Unknown
      QUADRUPOLE,                   ///< Quadrupole
      PAULIONTRAP,                  ///< Quadrupole ion trap / Paul ion trap
      RADIALEJECTIONLINEARIONTRAP,  ///< Radial ejection linear ion trap
      AXIALEJECTIONLINEARIONTRAP,   ///< Axial ejection linear ion trap
      TOF,                          ///< Time-of-flight
      SECTOR,                       ///< Magnetic sector
      FOURIERTRANSFORM,             ///< Fourier transform ion cyclotron resonance mass spectrometer
      IONSTORAGE,                   ///< Ion storage
      ESA,                          ///< Electrostatic energy analyzer
      IT,                           ///< Ion trap
      SWIFT,                        ///< Stored waveform inverse fourier transform
      CYCLOTRON,                    ///< Cyclotron
      ORBITRAP,                     ///< Orbitrap
      LIT,                          ///< Linear ion trap
      SIZE_OF_ANALYZERTYPE
    };
    /// Names of the analyzer types
    static const std::string NamesOfAnalyzerType[SIZE_OF_ANALYZERTYPE];

    /**
        @brief resolution method

        Which of the available standard measures is used to define whether two peaks are separate
    */
    enum ResolutionMethod
    {
      RESMETHNULL,                  ///< Unknown
      FWHM,                         ///< Full width at half max
      TENPERCENTVALLEY,             ///< Ten percent valley
      BASELINE,                     ///< Baseline
      SIZE_OF_RESOLUTIONMETHOD
    };
    /// Names of resolution methods
    static const std::string NamesOfResolutionMethod[SIZE_OF_RESOLUTIONMETHOD];

    /// Resolution type
    enum ResolutionType
    {
      RESTYPENULL,              ///< Unknown
      CONSTANT,                 ///< Constant
      PROPORTIONAL,             ///< Proportional
      SIZE_OF_RESOLUTIONTYPE
    };
    /// Names of resolution type
    static const std::string NamesOfResolutionType[SIZE_OF_RESOLUTIONTYPE];

    /// direction of scanning
    enum ScanDirection
    {
      SCANDIRNULL,              ///< Unknown
      UP,                       ///< Up
      DOWN,                     ///< Down
      SIZE_OF_SCANDIRECTION
    };
    /// Names of direction of scanning
    static const std::string NamesOfScanDirection[SIZE_OF_SCANDIRECTION];

    ///Scan law
    enum ScanLaw
    {
      SCANLAWNULL,              ///< Unknown
      EXPONENTIAL,              ///< Unknown
      LINEAR,                   ///< Linear
      QUADRATIC,                ///< Quadratic
      SIZE_OF_SCANLAW
    };
    /// Names of scan laws
    static const std::string NamesOfScanLaw[SIZE_OF_SCANLAW];

    ///Reflectron state
    enum ReflectronState
    {
      REFLSTATENULL,            ///< Unknown
      ON,                       ///< On
      OFF,                      ///< Off
      NONE,                     ///< None
      SIZE_OF_REFLECTRONSTATE
    };
    /// Names of reflectron states
    static const std::string NamesOfReflectronState[SIZE_OF_REFLECTRONSTATE];

    /// Constructor
    MassAnalyzer();
    /// Copy constructor
    MassAnalyzer(const MassAnalyzer &) = default;
    /// Move constructor
    MassAnalyzer(MassAnalyzer&&) = default;
    /// Destructor
    ~MassAnalyzer();

    /// Assignment operator
    MassAnalyzer & operator=(const MassAnalyzer&) = default;
    /// Move assignment operator
    MassAnalyzer& operator=(MassAnalyzer&&) & = default;

    /// Equality operator
    bool operator==(const MassAnalyzer & rhs) const;
    /// Equality operator
    bool operator!=(const MassAnalyzer & rhs) const;

    /// returns the analyzer type
    AnalyzerType getType() const;
    /// sets the analyzer type
    void setType(AnalyzerType type);

    /// returns the method used for determination of the resolution
    ResolutionMethod getResolutionMethod() const;
    /// sets the method used for determination of the resolution
    void setResolutionMethod(ResolutionMethod resolution_method);

    /// returns the resolution type
    ResolutionType getResolutionType() const;
    /// sets the resolution type
    void setResolutionType(ResolutionType resolution_type);

    /// returns the direction of scanning
    ScanDirection getScanDirection() const;
    /// sets the direction of scanning
    void setScanDirection(ScanDirection scan_direction);

    /// returns the scan law
    ScanLaw getScanLaw() const;
    /// sets the scan law
    void setScanLaw(ScanLaw scan_law);

    /// returns the reflectron state (for TOF)
    ReflectronState getReflectronState() const;
    /// sets the reflectron state (for TOF)
    void setReflectronState(ReflectronState reflecton_state);

    /**
        @brief returns the resolution

        The maximum m/z value at which two peaks can be resolved, according to one of the standard measures
    */
    double getResolution() const;
    /// sets the resolution
    void setResolution(double resolution);

    /// returns the mass accuracy i.e. how much the theoretical mass may differ from the measured mass (in ppm)
    double getAccuracy() const;
    /// sets the accuracy  i.e. how much the theoretical mass may differ from the measured mass  (in ppm)
    void setAccuracy(double accuracy);

    /// returns the scan rate (in s)
    double getScanRate() const;
    /// sets the scan rate (in s)
    void setScanRate(double scan_rate);

    /// returns the scan time for a single scan (in s)
    double getScanTime() const;
    /// sets the scan time for a single scan (in s)
    void setScanTime(double scan_time);

    /// returns the path length for a TOF mass analyzer (in meter)
    double getTOFTotalPathLength() const;
    /// sets the path length for a TOF mass analyzer (in meter)
    void setTOFTotalPathLength(double TOF_total_path_length);

    /// returns the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
    double getIsolationWidth() const;
    /// sets the isolation width i.e. in which m/z range the precursor ion is selected for MS to the n (in m/z)
    void setIsolationWidth(double isolation_width);

    /// returns the final MS exponent
    Int getFinalMSExponent() const;
    /// sets the final MS exponent
    void setFinalMSExponent(Int final_MS_exponent);

    /// returns the strength of the magnetic field (in T)
    double getMagneticFieldStrength() const;
    /// sets the strength of the magnetic field (in T)
    void setMagneticFieldStrength(double magnetic_field_strength);

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
    AnalyzerType type_;
    ResolutionMethod resolution_method_;
    ResolutionType resolution_type_;
    ScanDirection scan_direction_;
    ScanLaw scan_law_;
    ReflectronState reflectron_state_;
    double resolution_;
    double accuracy_;
    double scan_rate_;
    double scan_time_;
    double TOF_total_path_length_;
    double isolation_width_;
    Int final_MS_exponent_;
    double magnetic_field_strength_;
    Int order_;
  };
} // namespace OpenMS

