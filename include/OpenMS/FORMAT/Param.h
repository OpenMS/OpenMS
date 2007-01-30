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
		
		It is similar to a map<string,DataValue> but it also supports hierarchical 
		data and to save/load the contained data as XML.
		
		Hierachy levels are separated from each other and from the name by colons.
		
		Example: 'common:file_options:default_file_open_path = /share/'
		
		@note In the XML representation only the types 'int', 'string' ,'float' and 'double' are available.
		
		@see DefaultParamHandler
		
		@ingroup FileIO
	*/
	class Param
	{
	 public:
    /// Const iterator
    typedef std::map<std::string, DataValue>::const_iterator ConstIterator;

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

    /// Set a SignedInt value.
    void setValue(const std::string& key, SignedInt value);
    /// Set a float value.
    void setValue(const std::string& key, float value);
    /// Set a double value.
    void setValue(const std::string& key, double value);
    /// Set a string value.
    void setValue(const std::string& key, const std::string& value);

    /**
			 @brief Get a value by its key.

			 To check if there is no value for the given key, compare the return value with DataValue::EMPTY
    */
    const DataValue& getValue(const std::string& key) const;

    ///Returns the number of entries (leafs).
    UnsignedInt size() const;
    ///Returns if there are no entries.
    bool empty() const;
    /// Deletes all entries
    void clear();

    ///Insert all values of @p para and adds the prefix @p prefix.
    void insert(String prefix, const Param& para);
    ///Remove all entries that start with @p prefix.
    void remove(const std::string& prefix);

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
		Param copy(const std::string& prefix, bool remove_prefix=false, String new_prefix="") const;
	
		/** @brief Like copy(), but with support for "inherit" items.
				
				Inheritance is considered for "nodes" only, i.e. if old_prefix ends
				with ':'.  The old_prefix is <em>always</em> removed and replaced with
				new_prefix.  (Keeping old_prefix seems to make no sense in combination
				with inheritance.)
		 */
    Param copyWithInherit(const std::string& old_prefix, const std::string& new_prefix="") const;

    ///Write XML file.
    void store(const std::string& filename) const throw (Exception::UnableToCreateFile);
    ///Read XML file.
    void load(const std::string& filename) throw (Exception::FileNotFound,Exception::ParseError);

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
    void parseCommandLine(const int argc , char** argv, const std::map<std::string, std::string>& options_with_argument, const std::map<std::string, std::string>& options_without_argument, const std::string& misc="misc", const std::string& unknown="unknown");


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
    /// internal storage container
    std::map<std::string, DataValue> values_;

		friend std::ostream& operator << (std::ostream& os, const Param& param);

	 public:
		/**@brief Maximum number of inheritance steps allowed.
			 
		Usually you really won't care about this, thus I don't provide accessor functions. (Clemens)
		*/
		int inheritance_steps_max;

	};

	/**@brief Output of the object to a spream.
		 
	@relatesalso Param
	*/
	std::ostream& operator << (std::ostream& os, const Param& param);

} // namespace OpenMS

#endif // OPENMS_FORMAT_PARAM_H
