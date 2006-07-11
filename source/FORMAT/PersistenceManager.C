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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/PersistenceManager.h>
#include <OpenMS/FORMAT/RTTI.h>

#include <sstream>
#include <iostream>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{	
	PersistenceManager::PersistenceManager()
		: signature_constructor_(),
			todo_(),
			children_(),
			parents_(),
			current_(0)
	{
		//DPeak
		registerType_(RTTI::getStreamName< DPeak<1> >() , RTTI::getNew< DPeak<1> >);
		//DPickedPeak
		registerType_(RTTI::getStreamName< DPickedPeak<1> >() , RTTI::getNew< DPickedPeak<1> >);
		//DPeakArrayNonPolymorphic<1>
		registerType_(RTTI::getStreamName< DPeakArrayNonPolymorphic<1> >() , RTTI::getNew< DPeakArrayNonPolymorphic<1> >);
		//Spectrum<>
		registerType_(RTTI::getStreamName< MSSpectrum<> >() , RTTI::getNew< MSSpectrum<> >);
		//MSExperiment<>
		registerType_(RTTI::getStreamName< MSExperiment<> >() , RTTI::getNew< MSExperiment<> >);
	}
	
	PersistenceManager::~PersistenceManager()
	{
		
	}

	PersistenceManager& PersistenceManager::operator >>(PersistentObject*& object_ptr)
	{
		// if an error happened, just exit the loop to clean up the mess 
		bool error = false;
		// flag if this is the first object (it is the return value)
		bool first_object = true;
		// temporary pointer
		PersistentObject*	obj = 0;
		// streamName from RTTI
		String stream_name;


		// loop while there are objects to deserialize and no error occurred
		while (objectsToDeserialize() && !error) 
		{
			// retrieve the object signature
			getObjectHeader(stream_name);
				
			if (!signature_constructor_.has(stream_name)) 
			{
				error = true;
				break;
			} 
			
			// Create an instance of type_name 

			CreatePointerType	m = signature_constructor_[stream_name];
			obj = (PersistentObject*)m();

			// check whether the creation was successful
			if (obj == 0)
			{
				error = false;
				break;
			}
				
			// make the new object read itself
			obj->persistentRead(*this);

			// the first object is our return value
			if (first_object)
			{
				object_ptr = obj;
				first_object = false;
			}
		
			//deserialize subobjects
			while (!todo_.empty())
			{
				getObjectHeader(stream_name);
				//cout << endl <<"---- READING SUBOBJECT " << endl;
				const_cast< PersistentObject* >(todo_.begin()->first)->persistentRead(*this);
				
				//remove the subobject just read
				todo_.pop_front();
			}
		}

		if (error)
		{ 
			object_ptr = 0;
		}

		return *this;
	}

	void PersistenceManager::operator <<(PersistentObject& object)
	{
		//initialize with first object
		todo_.push_back(pair<PersistentObject*,string>(&object,"RootObject"));
		parents_.insert(pair<PersistentObject*,PersistentObject*>(&object,0));
		
		//do the actual main loop
		while (!todo_.empty())
		{
			todo_.front().first->persistentWrite(*this,todo_.front().second.c_str());
			todo_.pop_front();
		}

		//cleanup after work
		clear(); //cleanup of the derived class
		children_.clear();
		parents_.clear();
		current_ = 0;
	}
	
	/// ends writing the current object
	void 	PersistenceManager::writeObjectTrailer (const char* name )
	{
		writeTrailer(name);
	}

	void PersistenceManager::registerType_(const char* signature , const CreatePointerType create_pointer)
	{
		signature_constructor_.insert( make_pair( signature, create_pointer ) );
	}

}
