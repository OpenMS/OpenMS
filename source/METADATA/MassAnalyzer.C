// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MassAnalyzer.h>

using namespace std;

namespace OpenMS
{

	const std::string MassAnalyzer::NamesOfAnalyzerType[] = {"Unknown","Quadrupole","Quadrupole ion trap","Radial ejection linear ion trap","Axial ejection linear ion trap","Time-of-flight","Magnetic sector","Fourier transform ion cyclotron resonance mass spectrometer","Ion storage","Electrostatic energy analyzer","Ion trap","Stored waveform inverse fourier transform","Cyclotron","Orbitrap","Linear ion trap"};

	const std::string MassAnalyzer::NamesOfResolutionMethod[] = {"Unknown","Full width at half max","Ten percent valley","Baseline"};

	const std::string MassAnalyzer::NamesOfResolutionType[] = {"Unknown","Constant","Proportional"};

	const std::string MassAnalyzer::NamesOfScanDirection[] = {"Unknown","Up","Down"};

	const std::string MassAnalyzer::NamesOfScanLaw[] = {"Unknown","Exponential","Linar","Quadratic"};

	const std::string MassAnalyzer::NamesOfReflectronState[] = {"Unknown","On","Off","None"};
	
	MassAnalyzer::MassAnalyzer():
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
	
	MassAnalyzer::MassAnalyzer(const MassAnalyzer& source):
		MetaInfoInterface(source),
	  type_(source.type_),
	  resolution_method_(source.resolution_method_),
	  resolution_type_(source.resolution_type_),
	  scan_direction_(source.scan_direction_),
	  scan_law_(source.scan_law_),
	  reflectron_state_(source.reflectron_state_),
	  resolution_(source.resolution_),
	  accuracy_(source.accuracy_),
	  scan_rate_(source.scan_rate_),
	  scan_time_(source.scan_time_),
	  TOF_total_path_length_(source.TOF_total_path_length_),
	  isolation_width_(source.isolation_width_),
	  final_MS_exponent_(source.final_MS_exponent_),
	  magnetic_field_strength_(source.magnetic_field_strength_),
		order_(source.order_)
	{
	  
	}
	
	MassAnalyzer::~MassAnalyzer()
	{
	  
	}
	
	MassAnalyzer& MassAnalyzer::operator = (const MassAnalyzer& source)
	{
	  if (&source == this) return *this;
	  
	  order_ = source.order_;
    type_ = source.type_;
    resolution_method_ = source.resolution_method_;
    resolution_type_ = source.resolution_type_;
    scan_direction_ = source.scan_direction_;
    scan_law_ = source.scan_law_;
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
	 		order_ == rhs.order_ &&
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
	
	DoubleReal MassAnalyzer::getResolution() const 
	{
	  return resolution_; 
	}
	
	void MassAnalyzer::setResolution(DoubleReal resolution)
	{
	  resolution_ = resolution; 
	}
	
	DoubleReal MassAnalyzer::getAccuracy() const 
	{
	  return accuracy_; 
	}
	
	void MassAnalyzer::setAccuracy(DoubleReal accuracy)
	{
	  accuracy_ = accuracy; 
	}
	
	DoubleReal MassAnalyzer::getScanRate() const 
	{
	  return scan_rate_; 
	}
	
	void MassAnalyzer::setScanRate(DoubleReal scan_rate)
	{
	  scan_rate_ = scan_rate; 
	}
	
	DoubleReal MassAnalyzer::getScanTime() const 
	{
	  return scan_time_; 
	}
	
	void MassAnalyzer::setScanTime(DoubleReal scan_time)
	{
	  scan_time_ = scan_time; 
	}
	
	DoubleReal MassAnalyzer::getTOFTotalPathLength() const 
	{
	  return TOF_total_path_length_; 
	}
	
	void MassAnalyzer::setTOFTotalPathLength(DoubleReal TOF_total_path_length)
	{
	  TOF_total_path_length_ = TOF_total_path_length; 
	}
	
	DoubleReal MassAnalyzer::getIsolationWidth() const 
	{
	  return isolation_width_; 
	}
	
	void MassAnalyzer::setIsolationWidth(DoubleReal isolation_width)
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
	
	DoubleReal MassAnalyzer::getMagneticFieldStrength() const 
	{
	  return magnetic_field_strength_; 
	}
	
	void MassAnalyzer::setMagneticFieldStrength(DoubleReal magnetic_field_strength)
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

}

