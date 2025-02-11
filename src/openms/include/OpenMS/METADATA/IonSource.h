// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{

  /**
      @brief Description of an ion source (part of a MS Instrument)

      @ingroup Metadata
  */
  class OPENMS_DLLAPI IonSource :
    public MetaInfoInterface
  {
public:
    /// inlet type
    enum InletType
    {
      INLETNULL,                           ///< Unknown
      DIRECT,                              ///< Direct
      BATCH,                               ///< Batch (e.g. in MALDI)
      CHROMATOGRAPHY,                      ///< Chromatography (liquid)
      PARTICLEBEAM,                        ///< Particle beam
      MEMBRANESEPARATOR,                   ///< Membrane separator
      OPENSPLIT,                           ///< Open split
      JETSEPARATOR,                        ///< Jet separator
      SEPTUM,                              ///< Septum
      RESERVOIR,                           ///< Reservoir
      MOVINGBELT,                          ///< Moving belt
      MOVINGWIRE,                          ///< Moving wire
      FLOWINJECTIONANALYSIS,               ///< Flow injection analysis
      ELECTROSPRAYINLET,                   ///< Electro spray
      THERMOSPRAYINLET,                    ///< Thermo spray
      INFUSION,                            ///< Infusion
      CONTINUOUSFLOWFASTATOMBOMBARDMENT,   ///< Continuous flow fast atom bombardment
      INDUCTIVELYCOUPLEDPLASMA,            ///< Inductively coupled plasma
      MEMBRANE,                            ///< Membrane inlet
      NANOSPRAY,                           ///< Nanospray inlet
      SIZE_OF_INLETTYPE
    };
    /// Names of inlet types
    static const std::string NamesOfInletType[SIZE_OF_INLETTYPE];

    /// ionization method
    enum IonizationMethod
    {
      IONMETHODNULL,           ///< Unknown
      ESI,                     ///< electrospray ionisation
      EI,                      ///< electron ionization
      CI,                      ///< chemical ionisation
      FAB,                     ///< fast atom bombardment
      TSP,                     ///< thermospray
      LD,                      ///< laser desorption
      FD,                      ///< field desorption
      FI,                      ///< flame ionization
      PD,                      ///< plasma desorption
      SI,                      ///< secondary ion MS
      TI,                      ///< thermal ionization
      API,                     ///< atmospheric pressure ionisation
      ISI,                     ///<
      CID,                     ///< collision induced decomposition
      CAD,                     ///< collision activated decomposition
      HN,                      ///<
      APCI,                    ///< atmospheric pressure chemical ionization
      APPI,                    ///< atmospheric pressure photo ionization
      ICP,                     ///< inductively coupled plasma
      NESI,                    ///< Nano electrospray ionization
      MESI,                    ///< Micro electrospray ionization
      SELDI,                   ///< Surface enhanced laser desorption ionization
      SEND,                    ///< Surface enhanced neat desorption
      FIB,                     ///< Fast ion bombardment
      MALDI,                   ///< Matrix-assisted laser desorption ionization
      MPI,                     ///< Multiphoton ionization
      DI,                      ///< desorption ionization
      FA,                      ///< flowing afterglow
      FII,                     ///< field ionization
      GD_MS,                   ///< glow discharge ionization
      NICI,                    ///< negative ion chemical ionization
      NRMS,                    ///< neutralization reionization mass spectrometry
      PI,                      ///< photoionization
      PYMS,                    ///< pyrolysis mass spectrometry
      REMPI,                   ///< resonance enhanced multiphoton ionization
      AI,                      ///< adiabatic ionization
      ASI,                     ///< associative ionization
      AD,                      ///< autodetachment
      AUI,                     ///< autoionization
      CEI,                     ///< charge exchange ionization
      CHEMI,                   ///< chemi-ionization
      DISSI,                   ///< dissociative ionization
      LSI,                     ///< liquid secondary ionization
      PEI,                     ///< penning ionization
      SOI,                     ///< soft ionization
      SPI,                     ///< spark ionization
      SUI,                     ///< surface ionization
      VI,                      ///< vertical ionization
      AP_MALDI,                ///< atmospheric pressure matrix-assisted laser desorption ionization
      SILI,                    ///< desorption/ionization on silicon
      SALDI,                   ///< surface-assisted laser desorption ionization
      SIZE_OF_IONIZATIONMETHOD
    };
    /// Names of ionization methods
    static const std::string NamesOfIonizationMethod[SIZE_OF_IONIZATIONMETHOD];

    /// Polarity of the ion source
    enum Polarity
    {
      POLNULL,          ///< Unknown
      POSITIVE,         ///< Positive polarity
      NEGATIVE,         ///< Negative polarity
      SIZE_OF_POLARITY
    };
    /// Names of polarity of the ion source
    static const std::string NamesOfPolarity[SIZE_OF_POLARITY];

    /// Constructor
    IonSource();
    /// Copy constructor
    IonSource(const IonSource &) = default;
    /// Move constructor
    IonSource(IonSource&&) = default;
    /// Destructor
    ~IonSource();

    /// Assignment operator
    IonSource & operator=(const IonSource &) = default;
    /// Move assignment operator
    IonSource& operator=(IonSource&&) & = default;

    /// Equality operator
    bool operator==(const IonSource & rhs) const;
    /// Equality operator
    bool operator!=(const IonSource & rhs) const;

    /// returns the inlet type
    InletType getInletType() const;
    /// sets the  inlet type
    void setInletType(InletType inlet_type);

    /// returns the ionization method
    IonizationMethod getIonizationMethod() const;
    /// sets the ionization method
    void setIonizationMethod(IonizationMethod ionization_type);

    /// returns the ionization mode
    Polarity getPolarity() const;
    /// sets the ionization mode
    void setPolarity(Polarity polarity);

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
    InletType inlet_type_;
    IonizationMethod ionization_method_;
    Polarity polarity_;
    Int order_;

  };

} // namespace OpenMS

