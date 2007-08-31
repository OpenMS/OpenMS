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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_PARAM_H
#define OPENMS_FORMAT_PARAM_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <map>
#include <iostream>


namespace OpenMS
{     
	/**
		@brief Management and storage of INI files.
		
		This class provides a means to associate string names to int/double/string values.
		It also supports hierarchical data and to save/load the contained data as XML.
		Hierachy levels are separated from each other and from the name by colons. @n
		Example: 'common:file_options:default_file_open_path = /share/'
		
		In addition to the Type-Name-Value tuples descriptions can be added to each secection and value.
		See the setValue methods and setDescription(String). Newline characters in the description are
		possible.
		
		@note In the XML representation only the types 'int', 'string' ,'float' and 'double' are available.
		
		@see DefaultParamHandler
		
		@ingroup FileIO
	*/
	class Param
	{
	  public:
			/// Const iterator
			typedef std::map<String, DataValue>::const_iterator ConstIterator;
			
			/** @name Constructors and Destructors
			 */
			//@{
			/// Default construtor
			Param();
			
			/// Copy constructor
			Param(const Param& rhs);
			
			/// Destructor
			~Param();
			//@}
			
			/// Assignment operator
			Param& operator = (const Param& rhs);
			
			/// Equality operator
			bool operator == (const Param& rhs) const;
			
			/**
				@brief Set an Int value.
				
				@param key String key. Can contain ':' wich separated section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param user_parameter If 'true' this parameter is always shown. If 'false' it is only included in the expert mode
			*/
			void setValue(const String& key, Int value, const String& description="", bool user_parameter=false);
			
			/**
				@brief Set a float value.

				@param key String key. Can contain ':' wich separated section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param user_parameter If 'true' this parameter is always shown. If 'false' it is only included in the expert mode
			*/
			void setValue(const String& key, float value, const String& description="", bool user_parameter=false);
			
			/**
				@brief Set a double value.

				@param key String key. Can contain ':' wich separated section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param user_parameter If 'true' this parameter is always shown. If 'false' it is only included in the expert mode
			*/
			void setValue(const String& key, double value, const String& description="", bool user_parameter=false);
			
			/**
				@brief Set a string value.

				@param key String key. Can contain ':' wich separated section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param user_parameter If 'true' this parameter is always shown. If 'false' it is only included in the expert mode
			*/
			void setValue(const String& key, const String& value, const String& description="", bool user_parameter=false);
			
			/**
				@brief Get a value by its key.
				
				To check if there is no value for the given key, compare the return value with DataValue::EMPTY
			*/
			const DataValue& getValue(const String& key) const;
			
			/**
				@brief Sets a description for a key (section or actual value).
				
				@note The description is only set when a corresponding section or value exists.
			*/
			void setDescription(const String& location, const String& description);
			
			/**
				@brief Get a description by its key.
				
				To check if there is no description for the given key an empty string is returned.
			*/
			const String& getDescription(const String& key) const;
			
			///Returns the number of entries (leafs).
			UInt size() const;
			///Returns if there are no entries.
			bool empty() const;
			/// Deletes all entries
			void clear();
			
			///Insert all values of @p para and adds the prefix @p prefix.
			void insert(String prefix, const Param& para);
			///Remove all entries that start with @p prefix.
			void remove(const String& prefix);
			
			/**
				@brief Insert all values of @p para and adds the prefix @p prefix, if the values are not already set.
				
				@param defaults The default values. 
				@param prefix The prefix to add to all defaults. 
				@param showMessage If <tt>true</tt> each default that is actually set is printed to stdout as well.
				
				@see checkDefaults
			*/
			void setDefaults(const Param& defaults, String prefix="", bool showMessage=false);
			
			/**
				@brief Warns if a parameter is present for which no default value is specified.
				
				@param name A name that is displayed in error messages.
				@param defaults The default values. 
				@param prefix The prefix where to check for the defaults. 
				@param os The output stream for the warnings.
			*/
			void checkDefaults(const String& name, const Param& defaults, String prefix="", std::ostream& os = std::cout) const;
			
			/**
				 @brief Returns a new Param object containing all entries that start with @p prefix.
				
				 @param prefix should contain a ':' at the end if you want to extract a subtree.
				 Otherwise not only nodes, but as well values with that prefix are copied.
				 @param remove_prefix indicates if the prefix is removed before adding entries to the new Param
				 @param new_prefix is added to the front of all keys
			*/
			Param copy(const String& prefix, bool remove_prefix=false, String new_prefix="") const;
			
			/** @brief Like copy(), but with support for "inherit" items.
					
					Inheritance is considered for "nodes" only, i.e. if old_prefix ends
					with ':'.  The old_prefix is <em>always</em> removed and replaced with
					new_prefix.  (Keeping old_prefix seems to make no sense in combination
					with inheritance.)
			 */
			Param copyWithInherit(const String& old_prefix, const String& new_prefix="") const;
			
			///Write XML file.
			void store(const String& filename) const throw (Exception::UnableToCreateFile);
			///Read XML file.
			void load(const String& filename) throw (Exception::FileNotFound,Exception::ParseError);
			
			/**
				 @brief Parses command line arguments.
				
				 This method discriminates three types of arguments:<BR>
				 (1) options (starting with '-') that have a text argument<BR>
				 (2) options (starting with '-') that have no text argument<BR>
				 (3) text arguments (not starting with '-')
				
				 Command line arguments '-a avalue -b -c bvalue misc1 misc2' would be stored like this:<BR>
				 "prefix:-a" -> "avalue"<BR>
				 "prefix:-b" -> ""<BR>
				 "prefix:-c" -> "bvalue"<BR>
				 "prefix:misc" -> "misc1 misc2"<BR>
			
				 @param argc argc variable from command line
				 @param argv argv varaible from command line
				 @param prefix prefix for all options
			*/
			void parseCommandLine(const int argc , char** argv, String prefix = "");
			
			/**
				 @brief Parses command line arguments to specified key locations.
				
				 @param argc argc variable from command line
				 @param argv argv varaible from command line
				 @param options_with_argument a map of options that are followed by an argument (with key where they are stored)
				 @param options_without_argument a map of options that are not followed by an argument (with key where they are stored). Present options are set to the the string 'true'.
				 @param misc key where all non-option arguments are stored
				 @param unknown key where all unknown options are stored
			*/
			void parseCommandLine(const int argc , char** argv, const std::map<String, String>& options_with_argument, const std::map<String, String>& options_without_argument, const String& misc="misc", const String& unknown="unknown");
			
			
			/// Returns a constant iterator to the begin of the stored values.
			inline ConstIterator begin() const
			{
				return values_.begin();
			}
			
			/// Returns a constant iterator to the end of the stored values.
			inline ConstIterator end() const
			{
				return values_.end();
			}

		protected:
			/// internal storage containers
			std::map<String, DataValue> values_;
			std::map<String, String> descriptions_;
			
			friend std::ostream& operator << (std::ostream& os, const Param& param);

		public:
			/**
				@brief Maximum number of inheritance steps allowed.
					 
				Usually you really won't care about this, thus I don't provide accessor functions. (Clemens)
			*/
			Int inheritance_steps_max;

	};

	/**@brief Output of the object to a spream.
		 
	@relatesalso Param
	*/
	std::ostream& operator << (std::ostream& os, const Param& param);

} // namespace OpenMS

#endif // OPENMS_FORMAT_PARAM_H
