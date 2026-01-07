// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>

using namespace std;

namespace OpenMS
{

  const std::string MassAnalyzer::NamesOfAnalyzerType[] = {"Unknown", "Quadrupole", "Quadrupole ion trap", "Radial ejection linear ion trap", "Axial ejection linear ion trap", "Time-of-flight", "Magnetic sector", "Fourier transform ion cyclotron resonance mass spectrometer", "Ion storage", "Electrostatic energy analyzer", "Ion trap", "Stored waveform inverse fourier transform", "Cyclotron", "Orbitrap", "Linear ion trap"};

  const std::string MassAnalyzer::NamesOfResolutionMethod[] = {"Unknown", "Full width at half max", "Ten percent valley", "Baseline"};

  const std::string MassAnalyzer::NamesOfResolutionType[] = {"Unknown", "Constant", "Proportional"};

  const std::string MassAnalyzer::NamesOfScanDirection[] = {"Unknown", "Up", "Down"};

  const std::string MassAnalyzer::NamesOfScanLaw[] = {"Unknown", "Exponential", "Linar", "Quadratic"};

  const std::string MassAnalyzer::NamesOfReflectronState[] = {"Unknown", "On", "Off", "None"};

  MassAnalyzer::MassAnalyzer() :
    MetaInfoInterface(),
    type_(ANALYZERNULL),
    resolution_method_(RESMETHNULL),
    resolution_type_(RESTYPENULL),
    scan_direction_(SCANDIRNULL),
    scan_law_(SCANLAWNULL),
    reflectron_state_(REFLSTATENULL),
    resolution_(0.0),
    accuracy_(0.0),
    scan_rate_(0.0),
    scan_time_(0.0),
    TOF_total_path_length_(0.0),
    isolation_width_(0.0),
    final_MS_exponent_(0),
    magnetic_field_strength_(0.0),
    order_(0)
  {

  }

  MassAnalyzer::~MassAnalyzer() = default;

  bool MassAnalyzer::operator==(const MassAnalyzer & rhs) const
  {
    return order_ == rhs.order_ &&
           type_ == rhs.type_ &&
           resolution_method_ == rhs.resolution_method_ &&
           resolution_type_ == rhs.resolution_type_ &&
           scan_direction_ == rhs.scan_direction_ &&
           scan_law_ == rhs.scan_law_ &&
           reflectron_state_ == rhs.reflectron_state_ &&
           resolution_ == rhs.resolution_ &&
           accuracy_ == rhs.accuracy_ &&
           scan_rate_ == rhs.scan_rate_ &&
           scan_time_ == rhs.scan_time_ &&
           TOF_total_path_length_ == rhs.TOF_total_path_length_ &&
           isolation_width_ == rhs.isolation_width_ &&
           final_MS_exponent_ == rhs.final_MS_exponent_ &&
           magnetic_field_strength_ == rhs.magnetic_field_strength_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool MassAnalyzer::operator!=(const MassAnalyzer & rhs) const
  {
    return !(operator==(rhs));
  }

  MassAnalyzer::AnalyzerType MassAnalyzer::getType() const
  {
    return type_;
  }

  void MassAnalyzer::setType(MassAnalyzer::AnalyzerType type)
  {
    type_ = type;
  }

  MassAnalyzer::ResolutionMethod MassAnalyzer::getResolutionMethod() const
  {
    return resolution_method_;
  }

  void MassAnalyzer::setResolutionMethod(MassAnalyzer::ResolutionMethod resolution_method)
  {
    resolution_method_ = resolution_method;
  }

  MassAnalyzer::ResolutionType MassAnalyzer::getResolutionType() const
  {
    return resolution_type_;
  }

  void MassAnalyzer::setResolutionType(MassAnalyzer::ResolutionType resolution_type)
  {
    resolution_type_ = resolution_type;
  }

  MassAnalyzer::ScanDirection MassAnalyzer::getScanDirection() const
  {
    return scan_direction_;
  }

  void MassAnalyzer::setScanDirection(MassAnalyzer::ScanDirection scan_direction)
  {
    scan_direction_ = scan_direction;
  }

  MassAnalyzer::ScanLaw MassAnalyzer::getScanLaw() const
  {
    return scan_law_;
  }

  void MassAnalyzer::setScanLaw(MassAnalyzer::ScanLaw scan_law)
  {
    scan_law_ = scan_law;
  }

  MassAnalyzer::ReflectronState MassAnalyzer::getReflectronState() const
  {
    return reflectron_state_;
  }

  void MassAnalyzer::setReflectronState(MassAnalyzer::ReflectronState reflectron_state)
  {
    reflectron_state_ = reflectron_state;
  }

  double MassAnalyzer::getResolution() const
  {
    return resolution_;
  }

  void MassAnalyzer::setResolution(double resolution)
  {
    resolution_ = resolution;
  }

  double MassAnalyzer::getAccuracy() const
  {
    return accuracy_;
  }

  void MassAnalyzer::setAccuracy(double accuracy)
  {
    accuracy_ = accuracy;
  }

  double MassAnalyzer::getScanRate() const
  {
    return scan_rate_;
  }

  void MassAnalyzer::setScanRate(double scan_rate)
  {
    scan_rate_ = scan_rate;
  }

  double MassAnalyzer::getScanTime() const
  {
    return scan_time_;
  }

  void MassAnalyzer::setScanTime(double scan_time)
  {
    scan_time_ = scan_time;
  }

  double MassAnalyzer::getTOFTotalPathLength() const
  {
    return TOF_total_path_length_;
  }

  void MassAnalyzer::setTOFTotalPathLength(double TOF_total_path_length)
  {
    TOF_total_path_length_ = TOF_total_path_length;
  }

  double MassAnalyzer::getIsolationWidth() const
  {
    return isolation_width_;
  }

  void MassAnalyzer::setIsolationWidth(double isolation_width)
  {
    isolation_width_ = isolation_width;
  }

  Int MassAnalyzer::getFinalMSExponent() const
  {
    return final_MS_exponent_;
  }

  void MassAnalyzer::setFinalMSExponent(Int final_MS_exponent)
  {
    final_MS_exponent_ = final_MS_exponent;
  }

  double MassAnalyzer::getMagneticFieldStrength() const
  {
    return magnetic_field_strength_;
  }

  void MassAnalyzer::setMagneticFieldStrength(double magnetic_field_strength)
  {
    magnetic_field_strength_ = magnetic_field_strength;
  }

  Int MassAnalyzer::getOrder() const
  {
    return order_;
  }

  void MassAnalyzer::setOrder(Int order)
  {
    order_ = order;
  }

  StringList MassAnalyzer::getAllNamesOfAnalyzerType()
  {
    StringList names;
    names.reserve(SIZE_OF_ANALYZERTYPE);
    for (size_t i = 0; i < SIZE_OF_ANALYZERTYPE; ++i)
    {
      names.push_back(NamesOfAnalyzerType[i]);
    }
    return names;
  }

  StringList MassAnalyzer::getAllNamesOfResolutionMethod()
  {
    StringList names;
    names.reserve(SIZE_OF_RESOLUTIONMETHOD);
    for (size_t i = 0; i < SIZE_OF_RESOLUTIONMETHOD; ++i)
    {
      names.push_back(NamesOfResolutionMethod[i]);
    }
    return names;
  }

  StringList MassAnalyzer::getAllNamesOfResolutionType()
  {
    StringList names;
    names.reserve(SIZE_OF_RESOLUTIONTYPE);
    for (size_t i = 0; i < SIZE_OF_RESOLUTIONTYPE; ++i)
    {
      names.push_back(NamesOfResolutionType[i]);
    }
    return names;
  }

  StringList MassAnalyzer::getAllNamesOfScanDirection()
  {
    StringList names;
    names.reserve(SIZE_OF_SCANDIRECTION);
    for (size_t i = 0; i < SIZE_OF_SCANDIRECTION; ++i)
    {
      names.push_back(NamesOfScanDirection[i]);
    }
    return names;
  }

  StringList MassAnalyzer::getAllNamesOfScanLaw()
  {
    StringList names;
    names.reserve(SIZE_OF_SCANLAW);
    for (size_t i = 0; i < SIZE_OF_SCANLAW; ++i)
    {
      names.push_back(NamesOfScanLaw[i]);
    }
    return names;
  }

  StringList MassAnalyzer::getAllNamesOfReflectronState()
  {
    StringList names;
    names.reserve(SIZE_OF_REFLECTRONSTATE);
    for (size_t i = 0; i < SIZE_OF_REFLECTRONSTATE; ++i)
    {
      names.push_back(NamesOfReflectronState[i]);
    }
    return names;
  }

  const std::string& MassAnalyzer::analyzerTypeToString(AnalyzerType type)
  {
    if (type == SIZE_OF_ANALYZERTYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_ANALYZERTYPE");
    }
    return NamesOfAnalyzerType[static_cast<size_t>(type)];
  }

  MassAnalyzer::AnalyzerType MassAnalyzer::toAnalyzerType(const std::string& name)
  {
    auto first = &NamesOfAnalyzerType[0];
    auto last = &NamesOfAnalyzerType[SIZE_OF_ANALYZERTYPE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<AnalyzerType>(it - first);
  }

  const std::string& MassAnalyzer::resolutionMethodToString(ResolutionMethod method)
  {
    if (method == SIZE_OF_RESOLUTIONMETHOD)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_RESOLUTIONMETHOD");
    }
    return NamesOfResolutionMethod[static_cast<size_t>(method)];
  }

  MassAnalyzer::ResolutionMethod MassAnalyzer::toResolutionMethod(const std::string& name)
  {
    auto first = &NamesOfResolutionMethod[0];
    auto last = &NamesOfResolutionMethod[SIZE_OF_RESOLUTIONMETHOD];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<ResolutionMethod>(it - first);
  }

  const std::string& MassAnalyzer::resolutionTypeToString(ResolutionType type)
  {
    if (type == SIZE_OF_RESOLUTIONTYPE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_RESOLUTIONTYPE");
    }
    return NamesOfResolutionType[static_cast<size_t>(type)];
  }

  MassAnalyzer::ResolutionType MassAnalyzer::toResolutionType(const std::string& name)
  {
    auto first = &NamesOfResolutionType[0];
    auto last = &NamesOfResolutionType[SIZE_OF_RESOLUTIONTYPE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<ResolutionType>(it - first);
  }

  const std::string& MassAnalyzer::scanDirectionToString(ScanDirection direction)
  {
    if (direction == SIZE_OF_SCANDIRECTION)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_SCANDIRECTION");
    }
    return NamesOfScanDirection[static_cast<size_t>(direction)];
  }

  MassAnalyzer::ScanDirection MassAnalyzer::toScanDirection(const std::string& name)
  {
    auto first = &NamesOfScanDirection[0];
    auto last = &NamesOfScanDirection[SIZE_OF_SCANDIRECTION];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<ScanDirection>(it - first);
  }

  const std::string& MassAnalyzer::scanLawToString(ScanLaw law)
  {
    if (law == SIZE_OF_SCANLAW)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_SCANLAW");
    }
    return NamesOfScanLaw[static_cast<size_t>(law)];
  }

  MassAnalyzer::ScanLaw MassAnalyzer::toScanLaw(const std::string& name)
  {
    auto first = &NamesOfScanLaw[0];
    auto last = &NamesOfScanLaw[SIZE_OF_SCANLAW];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<ScanLaw>(it - first);
  }

  const std::string& MassAnalyzer::reflectronStateToString(ReflectronState state)
  {
    if (state == SIZE_OF_REFLECTRONSTATE)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value not allowed", "SIZE_OF_REFLECTRONSTATE");
    }
    return NamesOfReflectronState[static_cast<size_t>(state)];
  }

  MassAnalyzer::ReflectronState MassAnalyzer::toReflectronState(const std::string& name)
  {
    auto first = &NamesOfReflectronState[0];
    auto last = &NamesOfReflectronState[SIZE_OF_REFLECTRONSTATE];
    const auto it = std::find(first, last, name);
    if (it == last)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Value unknown", name);
    }
    return static_cast<ReflectronState>(it - first);
  }

}

