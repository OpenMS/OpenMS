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

#include <OpenMS/METADATA/HPLC.h>

using namespace std;

namespace OpenMS 
{
  HPLC::HPLC():
//  	PersistentObject(),
    instrument_(),
    column_(),
    temperature_(21),
    pressure_(0),
    flux_(0),
    comment_(),
    gradient_()
  {
    
  }
  
  HPLC::HPLC(const HPLC& source):
//    PersistentObject(source),
    instrument_(source.instrument_),
    column_(source.column_),
    temperature_(source.temperature_),
    pressure_(source.pressure_),
    flux_(source.flux_),
    comment_(source.comment_),
    gradient_(source.gradient_)
  {
    
  }
   
  HPLC::~HPLC()
  {
    
  }
  
  HPLC& HPLC::operator = (const HPLC& source)
  {
    if (source == *this) return *this;
    
//    PersistentObject::operator = (source);
    instrument_ = source.instrument_;
    column_ = source.column_;
    temperature_ = source.temperature_;
    pressure_ = source.pressure_;
    flux_ = source.flux_;
    comment_ = source.comment_;
    gradient_ = source.gradient_;
    
    return *this;
  }

  bool HPLC::operator == (const HPLC& rhs) const
  {
    return 
      ( instrument_ == rhs.instrument_ ) &&
      ( column_ == rhs.column_ ) &&
      ( temperature_ == rhs.temperature_ ) &&
      ( pressure_ == rhs.pressure_ ) &&
      ( flux_ == rhs.flux_ ) &&
      ( comment_ == rhs.comment_ ) &&
      ( gradient_ == rhs.gradient_ )    ;
  }

  bool HPLC::operator != (const HPLC& rhs) const
  {
    return !(operator == (rhs));
  }
 
  const String& HPLC::getInstrument() const
  {
    return instrument_;
  }
 
  void HPLC::setInstrument(const String& instrument)
  {
    instrument_ = instrument;
  }
 
 
  const String& HPLC::getColumn() const
  {
    return column_;
  }
 
  void HPLC::setColumn(const String& column)
  {
    column_ = column;
  }
 
 
  Int HPLC::getTemperature() const
  {
    return temperature_;
  }
 
  void HPLC::setTemperature(Int temperature)
  {
    temperature_ = temperature;
  }
 
 
  UInt HPLC::getPressure() const
  {
    return pressure_;
  }
 
  void HPLC::setPressure(UInt pressure)
  {
    pressure_ = pressure;
  }
 
 
  UInt HPLC::getFlux() const
  {
    return flux_;
  }
 
  void HPLC::setFlux(UInt flux)
  {
    flux_ = flux;
  }
 
 
  String HPLC::getComment() const
  {
    return comment_;
  }
 
  void HPLC::setComment(String comment)
  {
    comment_ = comment;
  }
 
 
  Gradient& HPLC::getGradient()
  {
    return gradient_;
  }
 
  const Gradient& HPLC::getGradient() const
  {
    return gradient_;
  }
 
  void HPLC::setGradient(const Gradient& gradient)
  {
    gradient_ = gradient;
  }
 
} // namespace OpenMS

