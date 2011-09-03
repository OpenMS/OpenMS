// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Anton Pervukhin <Anton.Pervukhin@CeBiTec.Uni-Bielefeld.DE> $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_EXCEPTION_H
#define OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_EXCEPTION_H

#include <string>
#include <stdexcept>

class ExceptionTest;

namespace OpenMS {

namespace ims {

/**
 * Base class for all ims-specific errors 
 */
class Exception : public std::exception {
public:
  explicit Exception() : msg("") { }
  explicit Exception(const std::string& message) : msg(message) { }
  virtual ~Exception() throw() { }

  /** @deprecated Use message() instead.
   * @todo think about deprecating this. the method is needed because
   * we derive from std::exception */
  virtual const char* what() const throw() {
    return this->msg.c_str();
  }

  virtual const std::string message() const {
    return this->msg;
  }
	
protected:
  void setMessage(const std::string& message) {
    this->msg = message;
  }
	
private:
  std::string msg;

	friend class ::ExceptionTest; // needs access to setMessage()
};

} // namespace ims

} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_MASSDECOMPOSITION_IMS_EXCEPTION_H
