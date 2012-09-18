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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_DB_PERSISTENTOBJECT_H
#define OPENMS_FORMAT_DB_PERSISTENTOBJECT_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>

namespace OpenMS
{ 
  /**
  	@brief Base class for all persistent objects.
  	
  	Interface for all classes that can be stored persistently in the OpenMS DB.

  	@ingroup DatabaseIO
  */
  class OPENMS_DLLAPI PersistentObject
  {

    public:
      /// Default constructor
      PersistentObject();

      /// Destructor
      virtual ~PersistentObject();
			
      /// Assignment operator
      PersistentObject& operator= (const PersistentObject& rhs);
			
      /**
      	@brief Returns the persistence id
      	
      	This id is only used in the DBAdapter the id is used to connect the object to the data stored in the DB.
      */
      const UID& getPersistenceId() const;

      /**
      	@brief Sets the persistence id
      	
      	This id is only used in the DBAdapter the id is used to connect the object to the data stored in the DB.
      	<BR>
      	Do not set the persistence id unless you know what you are doing!
      */
      void setPersistenceId(const UID& persistence_id);

      /**
      	@brief Clears the persistence id
      	
      	Sets the id to 0.<br>
      	@param deep determines which ids are cleared. <tt>false</tt> means that only the id of the current object is reset. 
      	<tt>true</tt> means that the ids of all sub-objects are reset as well (default).
      */
      void clearId(bool deep = true);

    protected:
    
      ///A persistence id used to refer the data back to the source
      UID persistence_id_;

      /**
      	@brief Clears the persistence id of all sub-objects.
      	
      	
      */
      virtual void clearChildIds_() =0;
  };

}
#endif
