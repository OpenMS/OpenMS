// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/CONCEPT/EnumHelpers.h>
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
    type_(AnalyzerType::ANALYZERNULL),
    resolution_method_(ResolutionMethod::RESMETHNULL),
    resolution_type_(ResolutionType::RESTYPENULL),
    scan_direction_(ScanDirection::SCANDIRNULL),
    scan_law_(ScanLaw::SCANLAWNULL),
    reflectron_state_(ReflectronState::REFLSTATENULL),
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
    return Helpers::enumNames(NamesOfAnalyzerType, AnalyzerType::SIZE_OF_ANALYZERTYPE);
  }

  StringList MassAnalyzer::getAllNamesOfResolutionMethod()
  {
    return Helpers::enumNames(NamesOfResolutionMethod, ResolutionMethod::SIZE_OF_RESOLUTIONMETHOD);
  }

  StringList MassAnalyzer::getAllNamesOfResolutionType()
  {
    return Helpers::enumNames(NamesOfResolutionType, ResolutionType::SIZE_OF_RESOLUTIONTYPE);
  }

  StringList MassAnalyzer::getAllNamesOfScanDirection()
  {
    return Helpers::enumNames(NamesOfScanDirection, ScanDirection::SIZE_OF_SCANDIRECTION);
  }

  StringList MassAnalyzer::getAllNamesOfScanLaw()
  {
    return Helpers::enumNames(NamesOfScanLaw, ScanLaw::SIZE_OF_SCANLAW);
  }

  StringList MassAnalyzer::getAllNamesOfReflectronState()
  {
    return Helpers::enumNames(NamesOfReflectronState, ReflectronState::SIZE_OF_REFLECTRONSTATE);
  }

  const std::string& MassAnalyzer::analyzerTypeToString(AnalyzerType type)
  {
    return Helpers::enumToName(NamesOfAnalyzerType, type, AnalyzerType::SIZE_OF_ANALYZERTYPE, "SIZE_OF_ANALYZERTYPE", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::AnalyzerType MassAnalyzer::toAnalyzerType(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfAnalyzerType, AnalyzerType::SIZE_OF_ANALYZERTYPE, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& MassAnalyzer::resolutionMethodToString(ResolutionMethod method)
  {
    return Helpers::enumToName(NamesOfResolutionMethod, method, ResolutionMethod::SIZE_OF_RESOLUTIONMETHOD, "SIZE_OF_RESOLUTIONMETHOD", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::ResolutionMethod MassAnalyzer::toResolutionMethod(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfResolutionMethod, ResolutionMethod::SIZE_OF_RESOLUTIONMETHOD, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& MassAnalyzer::resolutionTypeToString(ResolutionType type)
  {
    return Helpers::enumToName(NamesOfResolutionType, type, ResolutionType::SIZE_OF_RESOLUTIONTYPE, "SIZE_OF_RESOLUTIONTYPE", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::ResolutionType MassAnalyzer::toResolutionType(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfResolutionType, ResolutionType::SIZE_OF_RESOLUTIONTYPE, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& MassAnalyzer::scanDirectionToString(ScanDirection direction)
  {
    return Helpers::enumToName(NamesOfScanDirection, direction, ScanDirection::SIZE_OF_SCANDIRECTION, "SIZE_OF_SCANDIRECTION", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::ScanDirection MassAnalyzer::toScanDirection(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfScanDirection, ScanDirection::SIZE_OF_SCANDIRECTION, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& MassAnalyzer::scanLawToString(ScanLaw law)
  {
    return Helpers::enumToName(NamesOfScanLaw, law, ScanLaw::SIZE_OF_SCANLAW, "SIZE_OF_SCANLAW", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::ScanLaw MassAnalyzer::toScanLaw(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfScanLaw, ScanLaw::SIZE_OF_SCANLAW, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

  const std::string& MassAnalyzer::reflectronStateToString(ReflectronState state)
  {
    return Helpers::enumToName(NamesOfReflectronState, state, ReflectronState::SIZE_OF_REFLECTRONSTATE, "SIZE_OF_REFLECTRONSTATE", OPENMS_PRETTY_FUNCTION);
  }

  MassAnalyzer::ReflectronState MassAnalyzer::toReflectronState(const std::string& name)
  {
    return Helpers::enumFromName(NamesOfReflectronState, ReflectronState::SIZE_OF_REFLECTRONSTATE, name, "Value unknown", OPENMS_PRETTY_FUNCTION);
  }

}

