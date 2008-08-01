// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_DATASTRUCTURES_PARAM_H
#define OPENMS_DATASTRUCTURES_PARAM_H


#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <map>
#include <iostream>


namespace OpenMS
{	
	/**
		@brief Management and storage of parameters / INI files.
		
		This class provides a means to associate string names to int/double/string/StringList values.
		It allows for parameter hierarchies and to save/load the data as XML.
		Hierarchy levels are separated from each other by colons. @n
		Example: 'common:file_options:default_file_open_path = /share/'
		
		Each parameter and section has a description. Newline characters in the description are possible.
		
		Each parameter has an <i>advanced</i> flag that indicates if this parameter is shown to 
		all users (false) or in advanced mode only (true). This is mostly used in visualization.
		
		@todo Add support for string/int/float lists (Hiwi)
		
		@see DefaultParamHandler
		
		@ingroup Datastructures
	*/
	class Param
		: public Internal::XMLFile
	{
	  public:

			/// Parameter entry used to store the actual information inside of a Param entry
			struct ParamEntry
			{
				/// Default constructor
				ParamEntry();
				/// Constructor with name, description, value and advanced flag
				ParamEntry(const String& n, const DataValue& v, const String& d, bool a);
				/// Equality operator (only name and value are compared)
				bool operator==(const ParamEntry& rhs) const;
				
				/// Name of the entry
				String name;
				/// Description of the entry
				String description;
				/// Value associated with the entry
				DataValue value;
				/// Advanced parameter flag (If 'true' it is only shown in advanced mode)
				bool advanced;
				///@name Restrictions to accepted values (used in checkDefaults)
				//@{
				DoubleReal min_float; ///< Default: -std::numeric_limits<DoubleReal>::max()
				DoubleReal max_float; ///< Default: std::numeric_limits<DoubleReal>::max()
				Int min_int; ///< Default: -std::numeric_limits<Int>::max()
				Int max_int; ///< Default: std::numeric_limits<Int>::max()
				std::vector<String> valid_strings; ///< Default: empty
				//@}
			};
			
			///Node inside a Param object which is used to build the internal tree
			struct ParamNode
			{
				///Iterator for child nodes
				typedef std::vector<ParamNode>::iterator NodeIterator;
				///Iterator for entries
				typedef std::vector<ParamEntry>::iterator EntryIterator;
	
				///Default constructor
				ParamNode();						
				///Constructor with name and description
				ParamNode(const String& n, const String& d);
				///Equality operator (name, entries and subnodes are compared)
				bool operator==(const ParamNode& rhs) const;
				
				/**
					@brief Look up entry of this node (local search)
					
					Returns the end iterator if no entry is found
				*/
				EntryIterator findEntry(const String& name);
				/**
					@brief Look up subnode of this node (local search)
					
					Returns the end iterator if no entry is found
				*/
				NodeIterator findNode(const String& name);
				/**
					@brief Look up the parent node of the entry or node corresponding to @p name (tree search)
					
					Returns 0 if no entry is found
				*/
				ParamNode* findParentOf(const String& name);
				/**
					@brief Look up the entry corresponding to @p name (tree search)
					
					Returns 0 if no entry is found
				*/
				ParamEntry* findEntryRecursive(const String& name);
				
				///Inserts a @p node with the given @p prefix 
				void insert(const ParamNode& node, const String& prefix = "");
				///Inserts an @p entry with the given @p prefix 
				void insert(const ParamEntry& entry, const String& prefix = "");
				///Returns the number of entries in the whole subtree
				UInt size() const;
				///Returns the name suffix of a @p key (the part behind the last ':' character)
				String suffix(const String& key) const;
				
				/// Name of the node
				String name;
				/// Description of the node
				String description;
				/// Entries (leafs) in the node
				std::vector<ParamEntry> entries;
				/// Subnodes
				std::vector<ParamNode> nodes;
			};

		public:
		
			/// Forward const iterator for the Param class
			class ParamIterator
			{
				public:
					/// Struct that captures information on entered / left nodes for ParamIterator
					struct TraceInfo
					{
						inline TraceInfo(const String& n, const String& d, bool o)
							: name(n),
								description(d),
								opened(o)
						{	
						}
						/// name of the node
						String name;
						/// description of the node
						String description;
						/// If it was opened (true) or closed (false)
						bool opened;
					};

					///Default constructor used to create a past-the-end iterator
					ParamIterator();
					///Constructor for begin iterator
					ParamIterator(const Param::ParamNode& root);
					///Dereferencing
			    const Param::ParamEntry& operator*();
			    ///Dereferencing
		      const Param::ParamEntry* operator->();
			    ///Prefix increment operator 
			    ParamIterator& operator++();
			    ///Postfix increment operator 
			    ParamIterator operator++(Int);
		 			///Equality operator
			    bool operator==(const ParamIterator& rhs) const;
					///Equality operator
			    bool operator!=(const ParamIterator& rhs) const;
			    ///Returns the absolute path of the current element (including all sections)
			    String getName() const;
					///Returns the traceback of the opened and closed sections 
			    const std::vector< TraceInfo >& getTrace() const;
			    
				protected:
					///Pointer to the root node
					const Param::ParamNode* root_;
					///Index of the current ParamEntry (-1 means invalid)
					Int current_;
					///Pointers to the ParmNodes we are in
					std::vector<const Param::ParamNode*> stack_; 
					///Node traversal data during last ++ operation.
					std::vector< TraceInfo > trace_;
					
			};
			
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

			///@name Accessors for single parameters
			//@{
			/**
				@brief Set an Int value.
				
				@param key String key. Can contain ':' which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, Int value, const String& description="", bool advanced=false);
			/**
				@brief Set a UInt value.
				
				@param key String key. Can contain ':' which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, UInt value, const String& description="", bool advanced=false);
			/**
				@brief Set a float value.

				@param key String key. Can contain ':'which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, Real value, const String& description="", bool advanced=false);
			/**
				@brief Set a double value.

				@param key String key. Can contain ':' which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, DoubleReal value, const String& description="", bool advanced=false);
			/**
				@brief Set a string value.

				@param key String key. Can contain ':' which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, const String& value, const String& description="", bool advanced=false);
			/**
				@brief Set a StringList value.

				@param key String key. Can contain ':' which separates section names
				@param value The actual value
				@param description Verbose description of the parameter
				@param advanced If 'false' this parameter is always shown. If 'true' it is only shown in the expert mode
			*/
			void setValue(const String& key, const StringList& value, const String& description="", bool advanced=false);

			/**
				@brief Returns a value of a parameter.
			
				@exception Exception::ElementNotFound is thrown if the parameter does not exists.
			*/
			const DataValue& getValue(const String& key) const;
			/**
				@brief Returns the whole parameter entry.
			
				@exception Exception::ElementNotFound is thrown if the parameter does not exists.
			*/
			const ParamEntry& getEntry(const String& key) const;
			/**
				@brief Returns the description of a parameter.
			
				@exception Exception::ElementNotFound is thrown if the parameter does not exists.
			*/
			const String& getDescription(const String& key) const;
			/**
			  @brief Returns if the parameter is a advanced parameter.
			
				This is mainly used in the GUI to determine which parmeters are always displayed 
				and which parameters are displayed only in 'advanced mode'.
				
				@exception Exception::ElementNotFound is thrown if the parameter does not exists.
			*/
			bool isAdvancedParameter(const String& key) const;
			/**
				@brief Sets a description for an existing section
				
				Descriptions for values cannot be set with this method.
				They have to be set when inserting the value itself.
				
				@exception Exception::ElementNotFound is thrown if the section does not exists.
			*/
			void setSectionDescription(const String& key, const String& description);	
			/**
				@brief Returns the description corresponding to the section with name @p key.
			
				If the section does not exist an empty string is returned.
			*/
			const String& getSectionDescription(const String& key) const;
			/// Begin iterator for the internal tree
			ParamIterator begin() const;
			/// End iterator for the internal tree
			ParamIterator end() const;
			/// Tests if a parameter is set
			bool exists(const String& key) const;
			//@}
			
			///@name Manipulation of the whole parameter set
			//@{
			///Returns the number of entries (leafs).
			UInt size() const;
			///Returns if there are no entries.
			bool empty() const;
			/// Deletes all entries
			void clear();
			///Insert all values of @p param and adds the prefix @p prefix.
			void insert(String prefix, const Param& param);
			///Remove all entries that start with @p prefix.
			void remove(const String& prefix);
			/**
				@brief Returns a new Param object containing all entries that start with @p prefix.
				
				@param prefix should contain a ':' at the end if you want to extract a subtree.
							 Otherwise not only nodes, but as well values with that prefix are copied.
				@param remove_prefix indicates if the prefix is removed before adding entries to the new Param
			*/
			Param copy(const String& prefix, bool remove_prefix=false) const;
			//@}

			///@name Default value handling
			//@{
			/**
				@brief Insert all values of @p defaults and adds the prefix @p prefix, if the values are not already set.
				
				@param defaults The default values. 
				@param prefix The prefix to add to all defaults. 
				@param showMessage If <tt>true</tt> each default that is actually set is printed to stdout as well.
				
				@see checkDefaults
			*/
			void setDefaults(const Param& defaults, String prefix="", bool showMessage=false);
			/**
				@brief Checks the current parameter entries against given @p defaults
				
				Several checks are performed:
				- If a parameter is present for which no default value is specified, a warning is issued to @p os.
				- If the type of a parameter and its default do not match, an exception is thrown.
				- If a string parameter contains an invalid string, an exception is thrown.
				- If a numeric parameter is out of the valid range, an exception is thrown.
				
				@param name The name that is used in error messages.
				@param defaults The default values. 
				@param prefix The prefix where to check for the defaults. 
				@param os The output stream for the warnings.
				
				@exception Exception::InvalidParameter is thrown if errors occur during the check
			*/
			void checkDefaults(const String& name, const Param& defaults, String prefix="", std::ostream& os = std::cout) const;	
			/**
				@brief Sets the valid strings for the parameter @p key.
				
				It is only checked in checkDefaults(). 

				@exception Exception::InvalidParameter is thrown, if one of the strings contains a comma character
				@exception Exception::ElementNotFound exception is thrown, if the parameter is no string parameter
			*/
			void setValidStrings(const String& key, const std::vector<String>& strings);
			/**
				@brief Sets the minimum value for the integer parameter @p key. 
				
				It is only checked in checkDefaults(). 
				
				@exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
			*/			
			void setMinInt(const String& key, Int min);
			/**
				@brief Sets the maximum value for the integer parameter @p key. 
				
				It is only checked in checkDefaults().
				
				@exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
			*/
			void setMaxInt(const String& key, Int max);
			/**
				@brief Sets the minimum value for the floating point parameter @p key. 
				
				It is only checked in checkDefaults(). 
				
				@exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
			*/
			void setMinFloat(const String& key, DoubleReal min);
			/**
				@brief Sets the maximum value for the floating point parameter @p key. 
				
				It is only checked in checkDefaults(). 
				
				@exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
			*/
			void setMaxFloat(const String& key, DoubleReal max);
			//@}

			///@name Command line parsing
			//@{
			/**
				 @brief Parses command line arguments
				
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
			void parseCommandLine(const int argc , const char** argv, String prefix = "");
			/**
				 @brief Parses command line arguments to specified key locations.
				
				 @param argc argc variable from command line
				 @param argv argv varaible from command line
				 @param options_with_argument a map of options that are followed by an argument (with key where they are stored)
				 @param options_without_argument a map of options that are not followed by an argument (with key where they are stored). Present options are set to the the string 'true'.
				 @param misc key where all non-option arguments are stored (format: comma separated values enclosed in quotes, e.g. "val1","val2","val3","val with space")
				 @param unknown key where all unknown options are stored
			*/
			void parseCommandLine(const int argc , const char** argv, const std::map<String, String>& options_with_argument, const std::map<String, String>& options_without_argument, const String& misc="misc", const String& unknown="unknown");
			//@}
						
			///@name File I/O methods
			//@{
			/**
			  @brief Write XML file.

			  @exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename) const;
			/**
			  @brief Read XML file.

			  @exception Exception::FileNotFound is thrown if the file could not be found
			  @exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			void load(const String& filename);
			//@}
			
		protected:
			/**
			  @brief Returns a mutable reference to a parameter entry.

			  @exception Exception::ElementNotFound is thrown for unset parameters
			*/
			ParamEntry& getEntry_(const String& key) const;
	
			///Constructor from a node wich is used as root node
			Param(const Param::ParamNode& node);
			
			/// Invisible root node that stores all the data
			mutable Param::ParamNode root_;
	};

	///Output of Param to a stream.
	std::ostream& operator<< (std::ostream& os, const Param& param);


} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_PARAM_H
