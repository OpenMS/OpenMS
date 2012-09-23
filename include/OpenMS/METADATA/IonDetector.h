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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_IONDETECTOR_H
#define OPENMS_METADATA_IONDETECTOR_H

#include <OpenMS/METADATA/MetaInfoInterface.h>

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
    enum Type
    {
      TYPENULL,                                                                 ///< Unknown
      ELECTRONMULTIPLIER,                                           ///< Electron multiplier
      PHOTOMULTIPLIER,                                                  ///< Photo multiplier
      FOCALPLANEARRAY,                                                  ///< Focal plane array
      FARADAYCUP,                                                           ///< Faraday cup
      CONVERSIONDYNODEELECTRONMULTIPLIER,           ///< Conversion dynode electron multiplier
      CONVERSIONDYNODEPHOTOMULTIPLIER,                  ///< Conversion dynode photo multiplier
      MULTICOLLECTOR,                                                   ///< Multi-collector
      CHANNELELECTRONMULTIPLIER,                            ///< Channel electron multiplier
      CHANNELTRON,                                                          ///< channeltron
      DALYDETECTOR,                                                         ///< daly detector
      MICROCHANNELPLATEDETECTOR,                            ///< microchannel plate detector
      ARRAYDETECTOR,                                                    ///< array detector
      CONVERSIONDYNODE,                                                 ///< conversion dynode
      DYNODE,                                                                   ///< dynode
      FOCALPLANECOLLECTOR,                                          ///< focal plane collector
      IONTOPHOTONDETECTOR,                                          ///< ion-to-photon detector
      POINTCOLLECTOR,                                                   ///< point collector
      POSTACCELERATIONDETECTOR,                                 ///< postacceleration detector
      PHOTODIODEARRAYDETECTOR,                                  ///< photodiode array detector
      INDUCTIVEDETECTOR,                                            ///< inductive detector
      ELECTRONMULTIPLIERTUBE,                                   ///< electron multiplier tube
      SIZE_OF_TYPE
    };
    /// Names of detector types
    static const std::string NamesOfType[SIZE_OF_TYPE];

    /// Acquisition mode
    enum AcquisitionMode
    {
      ACQMODENULL,                          ///< Unknown
      PULSECOUNTING,                    ///< Pulse counting
      ADC,                                          ///< Analog-digital converter
      TDC,                                          ///< Time-digital converter
      TRANSIENTRECORDER,            ///< Transient recorder
      SIZE_OF_ACQUISITIONMODE
    };
    /// Names of acquisition modes
    static const std::string NamesOfAcquisitionMode[SIZE_OF_ACQUISITIONMODE];

    /// Constructor
    IonDetector();
    /// Copy constructor
    IonDetector(const IonDetector & source);
    /// Destructor
    ~IonDetector();

    /// Assignment operator
    IonDetector & operator=(const IonDetector & source);

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
    DoubleReal getResolution() const;
    /// sets the resolution (in ns)
    void setResolution(DoubleReal resolution);

    /// retruns the analog-to-digital converter sampling frequency (in Hz)
    DoubleReal getADCSamplingFrequency() const;
    /// sets the analog-to-digital converter sampling frequency (in Hz)
    void setADCSamplingFrequency(DoubleReal ADC_sampling_frequency);

    /**
        @brief returns the position of this part in the whole Instrument.

        Order can be ignored, as long the instrument has this default setup:
        - one ion source
        - one or many mass analyzers
        - one ion detector

        For more complex instuments, the order should be defined.
*/
    Int getOrder() const;
    /// sets the order
    void setOrder(Int order);

protected:
    Type type_;
    AcquisitionMode acquisition_mode_;
    DoubleReal resolution_;
    DoubleReal ADC_sampling_frequency_;
    Int order_;

  };
} // namespace OpenMS

#endif // OPENMS_METADATA_IONDETECTOR_H
