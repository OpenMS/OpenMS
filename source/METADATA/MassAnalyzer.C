// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MassAnalyzer.h>

using namespace std;

namespace OpenMS
{

	const std::string MassAnalyzer::NamesOfAnalyzerType[] = {"Unknown","QUADRUPOLE","PAULIONTRAP","RADIALEJECTIONLINEARIONTRAP","AXIALEJECTIONLINEARIONTRAP","TOF","SECTOR","FOURIERTRANSFORM","IONSTORAGE"};
	const std::string MassAnalyzer::NamesOfResolutionMethod[] = {"Unknown","FWHM","TENPERCENTVALLEY","BASELINE"};
	const std::string MassAnalyzer::NamesOfResolutionType[] = {"Unknown","CONSTANT","PROPORTIONAL"};
	const std::string MassAnalyzer::NamesOfScanFunction[] = {"Unknown","SELECTEDIONDETECTION","MASSSCAN"};
	const std::string MassAnalyzer::NamesOfScanDirection[] = {"Unknown","UP","DOWN"};
	const std::string MassAnalyzer::NamesOfScanLaw[] = {"Unknown","EXPONENTIAL","LINEAR","QUADRATIC"};
	const std::string MassAnalyzer::NamesOfTandemScanningMethod[] = {"Unknown","PRODUCTIONSCAN","PRECURSORIONSCAN","CONSTANTNEUTRALLOSS","SINGLEREACTIONMONITORING","MULTIPLEREACTIONMONITORING","SINGLEIONMONITORING","MULTIPLEIONMONITORING"};
	const std::string MassAnalyzer::NamesOfReflectronState[] = {"Unknown","ON","OFF","NONE"};
	
	MassAnalyzer::MassAnalyzer():
		MetaInfoInterface(),
		type_(ANALYZERNULL),
		resolution_method_(RESMETHNULL), 
		resolution_type_(RESTYPENULL),
		scan_function_(SCANFCTNULL),
		scan_direction_(SCANDIRNULL),
		scan_law_(SCANLAWNULL),
		tandem_scan_method_(TANDEMNULL),
		reflectron_state_(REFLSTATENULL),	
		resolution_(0.0),
		accuracy_(0.0),
		scan_rate_(0.0),
		scan_time_(0.0), 
		TOF_total_path_length_(0.0),
		isolation_width_(0.0),
	  final_MS_exponent_(0),
		magnetic_field_strength_(0.0)
	{
		
	}
	
	MassAnalyzer::MassAnalyzer(const MassAnalyzer& source):
		MetaInfoInterface(source),
	  type_(source.type_),
	  resolution_method_(source.resolution_method_),
	  resolution_type_(source.resolution_type_),
	  scan_function_(source.scan_function_),
	  scan_direction_(source.scan_direction_),
	  scan_law_(source.scan_law_),
	  tandem_scan_method_(source.tandem_scan_method_),
	  reflectron_state_(source.reflectron_state_),
	  resolution_(source.resolution_),
	  accuracy_(source.accuracy_),
	  scan_rate_(source.scan_rate_),
	  scan_time_(source.scan_time_),
	  TOF_total_path_length_(source.TOF_total_path_length_),
	  isolation_width_(source.isolation_width_),
	  final_MS_exponent_(source.final_MS_exponent_),
	  magnetic_field_strength_(source.magnetic_field_strength_)
	{
	  
	}
	
	MassAnalyzer::~MassAnalyzer()
	{
	  
	}
	
	MassAnalyzer& MassAnalyzer::operator = (const MassAnalyzer& source)
	{
	  if (&source == this) return *this;
	  
    type_ = source.type_;
    resolution_method_ = source.resolution_method_;
    resolution_type_ = source.resolution_type_;
    scan_function_ = source.scan_function_;
    scan_direction_ = source.scan_direction_;
    scan_law_ = source.scan_law_;
    tandem_scan_method_ = source.tandem_scan_method_;
    reflectron_state_ = source.reflectron_state_;
    resolution_ = source.resolution_;
    accuracy_ = source.accuracy_;
    scan_rate_ = source.scan_rate_;
    scan_time_ = source.scan_time_;
    TOF_total_path_length_ = source.TOF_total_path_length_;
    isolation_width_ = source.isolation_width_;
    final_MS_exponent_ = source.final_MS_exponent_;
    magnetic_field_strength_ = source.magnetic_field_strength_;
    MetaInfoInterface::operator=(source);
	    
	  return *this;
	}
	
  bool MassAnalyzer::operator== (const MassAnalyzer& rhs) const
  {
  	return 
	    type_ == rhs.type_ &&
	    resolution_method_ == rhs.resolution_method_ &&
	    resolution_type_ == rhs.resolution_type_ &&
	    scan_function_ == rhs.scan_function_ &&
	    scan_direction_ == rhs.scan_direction_ &&
	    scan_law_ == rhs.scan_law_ &&
	    tandem_scan_method_ == rhs.tandem_scan_method_ &&
	    reflectron_state_ == rhs.reflectron_state_ &&
	    resolution_ == rhs.resolution_ &&
	    accuracy_ == rhs.accuracy_ &&
	    scan_rate_ == rhs.scan_rate_ &&
	    scan_time_ == rhs.scan_time_ &&
	    TOF_total_path_length_ == rhs.TOF_total_path_length_ &&
	    isolation_width_ == rhs.isolation_width_ &&
	    final_MS_exponent_ == rhs.final_MS_exponent_ &&
	    magnetic_field_strength_ == rhs.magnetic_field_strength_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool MassAnalyzer::operator!= (const MassAnalyzer& rhs) const
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
	
	MassAnalyzer::ScanFunction MassAnalyzer::getScanFunction() const 
	{
	  return scan_function_; 
	}
	
	void MassAnalyzer::setScanFunction(MassAnalyzer::ScanFunction scan_function)
	{
	  scan_function_ = scan_function; 
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
	
	MassAnalyzer::TandemScanningMethod MassAnalyzer::getTandemScanMethod() const 
	{
	  return tandem_scan_method_; 
	}
	
	void MassAnalyzer::setTandemScanMethod(MassAnalyzer::TandemScanningMethod tandem_scan_method)
	{
	  tandem_scan_method_ = tandem_scan_method; 
	}
	
	MassAnalyzer::ReflectronState MassAnalyzer::getReflectronState() const 
	{
	  return reflectron_state_; 
	}
	
	void MassAnalyzer::setReflectronState(MassAnalyzer::ReflectronState reflectron_state)
	{
	  reflectron_state_ = reflectron_state; 
	}
	
	float MassAnalyzer::getResolution() const 
	{
	  return resolution_; 
	}
	
	void MassAnalyzer::setResolution(float resolution)
	{
	  resolution_ = resolution; 
	}
	
	float MassAnalyzer::getAccuracy() const 
	{
	  return accuracy_; 
	}
	
	void MassAnalyzer::setAccuracy(float accuracy)
	{
	  accuracy_ = accuracy; 
	}
	
	float MassAnalyzer::getScanRate() const 
	{
	  return scan_rate_; 
	}
	
	void MassAnalyzer::setScanRate(float scan_rate)
	{
	  scan_rate_ = scan_rate; 
	}
	
	float MassAnalyzer::getScanTime() const 
	{
	  return scan_time_; 
	}
	
	void MassAnalyzer::setScanTime(float scan_time)
	{
	  scan_time_ = scan_time; 
	}
	
	float MassAnalyzer::getTOFTotalPathLength() const 
	{
	  return TOF_total_path_length_; 
	}
	
	void MassAnalyzer::setTOFTotalPathLength(float TOF_total_path_length)
	{
	  TOF_total_path_length_ = TOF_total_path_length; 
	}
	
	float MassAnalyzer::getIsolationWidth() const 
	{
	  return isolation_width_; 
	}
	
	void MassAnalyzer::setIsolationWidth(float isolation_width)
	{
	  isolation_width_ = isolation_width; 
	}
	
	SignedInt MassAnalyzer::getFinalMSExponent() const 
	{
	  return final_MS_exponent_; 
	}
	
	void MassAnalyzer::setFinalMSExponent(SignedInt final_MS_exponent)
	{
	  final_MS_exponent_ = final_MS_exponent; 
	}
	
	float MassAnalyzer::getMagneticFieldStrength() const 
	{
	  return magnetic_field_strength_; 
	}
	
	void MassAnalyzer::setMagneticFieldStrength(float magnetic_field_strength)
	{
	  magnetic_field_strength_ = magnetic_field_strength; 
	}
	
}

