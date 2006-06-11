// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: PersistenceManager.h,v 1.18 2006/03/14 16:59:37 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PERSISTENCEMANAGER_H
#define OPENMS_FORMAT_PERSISTENCEMANAGER_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/FORMAT/RTTI.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>

#include <string>
#include <list>
#include <map>
#include <vector>
#include <typeinfo>
#include <iostream>

namespace OpenMS
{
  class PersistentObject;
  
	/**
		@brief Base class for persistent storage of objects.
		
		All classes that implement the PersistentObject interface can be serialized and deserialized.
		
		@note Reconstructing complicated pointer structures is not yet supported. Only references work so far.
		
  	@ingroup Format
  	
  	@todo write test <BR> add FloatKernelTraits <BR> add Header/Trailer for object arrays (needed e.g. for XML)
	*/
  class PersistenceManager
  {
    public:
    	/// Pointer to a create method
    	typedef void* (*CreatePointerType)();

    	/// Default constructor
      PersistenceManager();
      /// Destructor
      virtual ~PersistenceManager();
			
			/**
				@brief Writes a persistent object. This is the only method a user should call for writing.
				
				In some implementations e.g. DBAdapter the persistent UID of @p object is set to identify the object.				
			*/
			void operator <<(PersistentObject& object);

			/**	
				@brief Reads a persistent object. This is the only method a user should call for reading.
				
				In some derived classes e.g. DBAdapter it is unclear what object should be read.
				In that case the derived class provieds a method to identify the object by some identifier.
				If no identifier was provided, @p object_ptr will be 0.
			*/
			PersistenceManager& operator >> (PersistentObject*& object_ptr);
      
			/** 
				@name Layer 1 Methods for writing
				
				Methods which are used inside the PersistentObject classes for writing
			*/
			//@{
			/// Begins writing the current object
			template<typename T> void writeObjectHeader (const T* object, const char* name=0)
      {
      	current_ = object;
      	writeHeader(streamClassName(typeid(*object)).c_str(), name, object);
      }
			
			/// Write a primitive data type i.e. SignedInt, UnsingedInt, double, string.
      template<typename T> void writePrimitive(const T& value, const char* name)
      {
      	writePrimitiveHeader(streamClassName(typeid(value)).c_str(), name);
				put(value);
				writePrimitiveTrailer();
      }
      
      /// Write an object reference 
      template<typename T> void writeObjectReference (const T& object, const char* name)
      {
      	const PersistentObject* tmp((const PersistentObject*) &object);
      	//std::cout << "writeObjectReference: " << tmp << " name: " << name << " parent: " << current_ << std::endl;
				if (&object != 0)
				{
					todo_.push_back( std::make_pair(tmp,name) );
					children_.insert(std::make_pair(current_,tmp));
					parents_.insert(std::make_pair(tmp,current_));
				}				
      }
      
      /// Write an array of object references
      template<typename T> void writeObjectArray (const T& array, const char* name, Size size) 
      {
      	//std::cout << "writeObjectArray: "<< PointerSizeInt(&array) <<" (size: " << size << ")" <<std::endl;
				for (Position i = 0; i < size; i++)
				{
					writeObjectReference(array[i],name);
				}
	    }
			
			/// Ends writing the current object
			void writeObjectTrailer (const char* name);
			

			/** 
				@name Layer 1 Methods for reading
				
				These methods used inside the PersistentObject classes for reading
			*/
			//@{	

			/**	
				@brief Reads a primitive data type e.g. int, float, string.
				
				A mutable refernence to the primitive and the its name is given and the primitive is filled with
				the associated value.
				
				@return	true if reading was successful
			*/
			template <typename T> bool readPrimitive(T& t, const char* name)
			{
				if (!checkPrimitiveHeader(RTTI::getStreamName<T>(), name))
				{
					return false;
				}
				get(t);
				return checkPrimitiveTrailer();
			}

			template <typename T>
			bool readObjectReference(T& object, const char* name)
			{
				if (!checkObjectReferenceHeader(RTTI::getStreamName<T>(), name))
				{
					return false;
				}

				todo_.push_back( std::make_pair( (const PersistentObject*) &object , name ) );
		
				return checkPrimitiveTrailer();
			}
			
			//@}	

			/** 
				@name Layer 0 Methods
				
				These methods are all purely virtual and must be implemented in the derived classes
			*/
			//@{
			
			///Header that is called at the beginning of each object
			virtual void writeHeader (const char* signature, const char* name, const PersistentObject* object ) =0;
			
			///Trailer that is called at the end of each object
			virtual void writeTrailer (const char* name) =0;
			
			///Header that is called at the end of each primitive
      virtual void writePrimitiveHeader (const char* signature, const char* name) =0;
			
			///Trailer that is called at the end of each primitive
			virtual void writePrimitiveTrailer() =0;
			
			///Writes a signed integer primitive
			virtual void put(const SignedInt value) =0;
			
			///Writes an unsigned integer primitive
			virtual void put(const UnsignedInt value) =0;
			
			///Writes a floating point primitive
			virtual void put(const double value) =0;
			
			///Writes a string primitive
			virtual void put(const std::string& value) =0;

			///Do the cleanup after all is done
			virtual void clear() =0;

			///Get an (unknown) object header. The stream name of the object is returned in @p stream_name.
			virtual bool getObjectHeader(String& stream_name) = 0;
			
			///Returns true if there are objects to serialize left
			virtual bool objectsToDeserialize() =0;

			/**	
				@brief Check for a type header and name for a primitive type.
				
				@param	type_name the stream name of the primitive
				@param	name the name of the primitive
				@return	true if type and name of the primitive match
			*/
			virtual bool checkPrimitiveHeader(const char* stream_name, const char* name) = 0;

			/**	
				@brief Check for header for a reference to a PersistentObject.
				
				@param	type_name the stream name of the object type
				@param	name the name of the object
				@return	true if the header was correct
			*/
			virtual bool checkObjectReferenceHeader(const char* type_name, const char* name) = 0;

			/**	
				@brief Check for the trailer of a primitive type.
				
				
				@return true if the trailer was correct
			*/
			virtual bool checkPrimitiveTrailer() = 0;

			///Read a double from the input stream.
			virtual void get(double& d) = 0;

			///Read a SignedInt from the input stream.
			virtual void get(UnsignedInt& i) = 0;
	
			///Read a SignedInt from the input stream.
			virtual void get(SignedInt& i) = 0;
	
			///Read a string from the output stream.
			virtual void get(string& s) = 0;

			///Read a UID from the output stream.
			virtual void get(UID& id) = 0;
			//@}

    protected:
    	///All types that need to be deserialized have to be registered with this method in the constructor.
    	void registerType_(const char* signature , const CreatePointerType create_pointer);
			
			/// Mapping of the RTTI stream name of an object to a creator method.
			HashMap<std::string, CreatePointerType> signature_constructor_;
			
      /// list of objects to process (with names)
      std::list< std::pair<const PersistentObject*, std::string> > todo_;

      ///a map that stored the child connections between objects
      std::multimap<const PersistentObject*, const PersistentObject* > children_;
      
      ///a map that stored the parent connections between objects
      std::map<const PersistentObject*, const PersistentObject* > parents_;
      	
      ///Pointer to the current object
      const PersistentObject* current_;
  };

}
#endif
