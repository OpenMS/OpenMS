// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <set>

namespace OpenMS
{

  template <class Key, class T>
  class Map;

  namespace Logger
  {
    class LogStream;
  }

  /**
    @brief Management and storage of parameters / INI files.

    This class provides a means to associate string names to int/double/string/StringList values.
    It allows for parameter hierarchies and to save/load the data as XML.
    Hierarchy levels are separated from each other by colons. @n
    Example: 'common:file_options:default_file_open_path = /share/'

    Each parameter and section has a description. Newline characters in the description are possible.

    Each parameter can be annotated with an arbitrary number of tags. Tags must not contain comma characters!
    @n E.g. the <i>advanced</i> tag indicates if this parameter is shown to all users or in advanced mode only.

    @see DefaultParamHandler

    @ingroup Datastructures
  */
  class OPENMS_DLLAPI Param
  {
public:

    /// Parameter entry used to store the actual information inside of a Param entry
    struct OPENMS_DLLAPI ParamEntry
    {
      /// Default constructor
      ParamEntry();
      /// Constructor with name, description, value and advanced flag
      ParamEntry(const String& n, const DataValue& v, const String& d, const StringList& t = StringList());
      /// Copy constructor
      ParamEntry(const ParamEntry&) = default;
      /// Move constructor
      ParamEntry(ParamEntry&&) = default;
      /// Destructor
      ~ParamEntry();

      /// Assignment operator
      ParamEntry& operator=(const ParamEntry&) = default;
      /// Move assignment operator
      ParamEntry& operator=(ParamEntry&&) & = default;

      /// Check if 'value' fulfills restrictions
      bool isValid(String& message) const;
      /// Equality operator (only name and value are compared)
      bool operator==(const ParamEntry& rhs) const;

      /// Name of the entry
      String name;
      /// Description of the entry
      String description;
      /// Value associated with the entry
      DataValue value;
      /// Tags list, used e.g. for advanced parameter tag
      std::set<String> tags;
      ///@name Restrictions to accepted values (used in checkDefaults)
      //@{
      double min_float; ///< Default: - std::numeric_limits<double>::max()
      double max_float; ///< Default: std::numeric_limits<double>::max()
      Int min_int; ///< Default: - std::numeric_limits<Int>::max()
      Int max_int; ///< Default: std::numeric_limits<Int>::max()
      std::vector<String> valid_strings; ///< Default: empty
      //@}
    };

    ///Node inside a Param object which is used to build the internal tree
    struct OPENMS_DLLAPI ParamNode
    {
      ///Iterator for child nodes
      typedef std::vector<ParamNode>::iterator NodeIterator;
      ///Iterator for entries
      typedef std::vector<ParamEntry>::iterator EntryIterator;
      ///Iterator for child nodes
      typedef std::vector<ParamNode>::const_iterator ConstNodeIterator;
      ///Iterator for entries
      typedef std::vector<ParamEntry>::const_iterator ConstEntryIterator;

      /// Default constructor
      ParamNode();
      /// Constructor with name and description
      ParamNode(const String& n, const String& d);
      /// Copy constructor
      ParamNode(const ParamNode&) = default;
      /// Move constructor
      ParamNode(ParamNode&&) = default;
      /// Destructor
      ~ParamNode();

      /// Assignment operator
      ParamNode& operator=(const ParamNode&) = default;
      /// Move assignment operator
      ParamNode& operator=(ParamNode&&) & = default;

      /// Equality operator (name, entries and subnodes are compared)
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
      Size size() const;
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
    class OPENMS_DLLAPI ParamIterator
    {
public:
      /// Struct that captures information on entered / left nodes for ParamIterator
      struct OPENMS_DLLAPI TraceInfo
      {
        /// Constructor with name, description, and open flag
        inline TraceInfo(const String& n, const String& d, bool o) :
          name(n),
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

      /// Default constructor used to create a past-the-end iterator
      ParamIterator();
      /// Constructor for begin iterator
      ParamIterator(const Param::ParamNode& root);
      /// Destructor
      ~ParamIterator();
      /// Dereferencing
      const Param::ParamEntry& operator*();
      /// Dereferencing
      const Param::ParamEntry* operator->();
      /// Prefix increment operator
      ParamIterator& operator++();
      /// Postfix increment operator
      ParamIterator operator++(int);
      /// Equality operator
      bool operator==(const ParamIterator& rhs) const;
      /// Equality operator
      bool operator!=(const ParamIterator& rhs) const;
      /// Returns the absolute path of the current element (including all sections)
      String getName() const;
      /// Returns the traceback of the opened and closed sections
      const std::vector<TraceInfo>& getTrace() const;

protected:
      /// Pointer to the root node
      const Param::ParamNode* root_;
      /// Index of the current ParamEntry (-1 means invalid)
      Int current_;
      /// Pointers to the ParamNodes we are in
      std::vector<const Param::ParamNode*> stack_;
      /// Node traversal data during last ++ operation.
      std::vector<TraceInfo> trace_;

    };

    /// Default constructor
    Param();

    /// Copy constructor
    Param(const Param&) = default;

    /// Move constructor
    Param(Param&&) = default;

    /// Destructor
    ~Param();

    /// Assignment operator
    Param& operator=(const Param&) = default;

    /// Move assignment operator
    Param& operator=(Param&&) & = default;

    /// Equality operator
    bool operator==(const Param& rhs) const;

    /// Begin iterator for the internal tree
    ParamIterator begin() const;

    /// End iterator for the internal tree
    ParamIterator end() const;

    ///@name Accessors for single parameters
    //@{

    /**
      @brief Sets a value.

      @param key String key. Can contain ':' which separates section names
      @param value The actual value
      @param description Verbose description of the parameter
      @param tags list of tags associated to this parameter
    */
    void setValue(const String& key, const DataValue& value, const String& description = "", const StringList& tags = StringList());

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
      @brief Tests if a parameter is set (expecting its fully qualified name, e.g., TextExporter:1:proteins_only)

      @param key The fully qualified name of the parameter to check.
      @return True if the parameter exists, false otherwise.
    */
    bool exists(const String& key) const;

    /**
      @brief Find leaf node by name (if it exists).

      @param leaf The name of the parameter to find excluding the path parameter, e.g., given the parameter TextExporter:1:proteins_only the leaf would be named proteins_only.
      @return Returns end() if leaf does not exist.
    */
    ParamIterator findFirst(const String& leaf) const;

    /**
      @brief Find next leaf node by name (if it exists), not considering the @p start_leaf

      @param leaf The name of the parameter to find excluding the path parameter, e.g., given the parameter TextExporter:1:proteins_only the leaf would be named proteins_only.
      @param start_leaf The already found leaf, that should not be considered during this search.
      @return Returns end() if leaf does not exist.
    */
    ParamIterator findNext(const String& leaf, const ParamIterator& start_leaf) const;
    //@}

    ///@name Tags handling
    //@{

    /**
      @brief Adds the tag @p tag to the entry @p key

      E.g. "advanced", "required", "input file", "output file"

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
      @exception Exception::InvalidValue is thrown if the tag contain a comma character.
    */
    void addTag(const String& key, const String& tag);

    /**
      @brief Adds the tags in the list @p tags to the entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
      @exception Exception::InvalidValue is thrown if a tag contain a comma character.
    */
    void addTags(const String& key, const StringList& tags);

    /**
      @brief Returns if the parameter @p key has a tag

      Example: The tag 'advanced' is used in the GUI to determine which parameters are always displayed
      and which parameters are displayed only in 'advanced mode'.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    bool hasTag(const String& key, const String& tag) const;

    /**
      @brief Returns the tags of entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    StringList getTags(const String& key) const;

    /**
      @brief Removes all tags from the entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    void clearTags(const String& key);
    //@}


    ///@name Descriptions handling
    //@{

    /**
      @brief Returns the description of a parameter.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    const String& getDescription(const String& key) const;

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

    /**
    @brief Adds a parameter section under the path @p key with the given @p description.

    If the section already exists, the description is only overwritten if not empty.
    */
    void addSection(const String& key, const String& description);
    //@}

    ///@name Manipulation of the whole parameter set
    //@{

    ///Returns the number of entries (leafs).
    Size size() const;

    ///Returns if there are no entries.
    bool empty() const;

    /// Deletes all entries
    void clear();

    /// Insert all values of @p param and adds the prefix @p prefix.
    /// You should append ':' to prefix manually when you want it to be a section.
    void insert(const String& prefix, const Param& param);

    /**
      @brief Remove the entry @p key or a section @p key (when suffix is ':')

      Remove deletes either an entry or a section (when @p key ends with ':'),
      by matching the exact name. No partial matches are accepted.

      If an empty internal node remains, the tree is pruned until every node has either a successor node
      or a leaf, i.e. no naked nodes remain.
    */
    void remove(const String& key);

    /**
      @brief Remove all entries that start with @p prefix

      Partial are valid as well. All entries and sections which match the prefix are deleted.

      If an empty internal node remains, the tree is pruned until every node has either a successor node
      or a leaf, i.e. no naked nodes remain.
    */
    void removeAll(const String& prefix);

    /**
      @brief Returns a new Param object containing all entries that start with @p prefix.

      @param prefix should contain a ':' at the end if you want to extract a subtree.
             Otherwise not only nodes, but as well values with that prefix are copied.
      @param remove_prefix indicates if the prefix is removed before adding entries to the new Param
    */
    Param copy(const String& prefix, bool remove_prefix = false) const;

    /**
      @brief Returns a new Param object containing all entries in the given subset.

      @param subset The subset of Param nodes that should be copied out of the object
             here. Includes values etc. Does not check any compatibility. Just matches the names.
      @note Only matches entries and nodes at the root=top level and copies over whole subtrees if matched.
            This function is mainly used for copying subsection parameters that were not registered as
            as an actual subsection e.g. for backwards compatibility of param names.
    */
    Param copySubset(const Param& subset) const;

    /**
      @brief Rescue parameter <b>values</b> from @p p_outdated to current param

      Calls ::update(p_outdated, true, add_unknown, false, false, OPENMS_LOG_WARN) and returns its value.
    */
    bool update(const Param& p_outdated, const bool add_unknown = false);

    /**
      @brief Rescue parameter <b>values</b> from @p p_outdated to current param

      Calls ::update(p_outdated, true, add_unknown, false, false, stream) and returns its value.
    */
    bool update(const Param& p_outdated, const bool add_unknown, Logger::LogStream& stream);
   
   
    /**
      @brief Rescue parameter <b>values</b> from @p p_outdated to current param

      All parameters present in both param objects will be transferred into this object, given that:
      <ul>
        <li>the name is equal</li>
        <li>the type is equal</li>
        <li>the value from @p p_outdated meets the new restrictions</li>
      </ul>

      Not transferred are values from parameter "version" (to preserve the new version) or "type" (to preserve layout).

      @param p_outdated Old/outdated param object, whose values (as long as they are still valid) are used to update this object 
      @param verbose Print information about expected value updates
      @param add_unknown Add unknown parameters from @p_outdated to this param object.
      @param fail_on_invalid_values Return false if outdated parameters hold invalid values
      @param fail_on_unknown_parameters Return false if outdated parameters contain unknown parameters (takes precedence over @p add_unknown)
      @param stream The stream where all the logging output is send to.
      @return true on success, false on failure
    */
    bool update(const Param& p_outdated, bool verbose, bool add_unknown, bool fail_on_invalid_values, bool fail_on_unknown_parameters, Logger::LogStream& stream);

    /**
      @brief Adds missing parameters from the given param @p toMerge to this param. Existing parameters will not be modified.

      @param toMerge The Param object from which parameters should be added to this param.
    */
    void merge(const Param& toMerge);

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
    void setDefaults(const Param& defaults, const String& prefix = "", bool showMessage = false);

    /**
      @brief Checks the current parameter entries against given @p defaults

      Several checks are performed:
      - If a parameter is present for which no default value is specified, a warning is issued to @p os.
      - If the type of a parameter and its default do not match, an exception is thrown.
      - If a string parameter contains an invalid string, an exception is thrown.
      -	If parameter entry is a string list, an exception is thrown, if one or more list members are invalid strings
      - If a numeric parameter is out of the valid range, an exception is thrown.
      - If entry is a numeric list an exception is thrown, if one or more list members are out of the valid range

      @param name The name that is used in error messages.
      @param defaults The default values.
      @param prefix The prefix where to check for the defaults.

      Warnings etc. will be send to OPENMS_LOG_WARN.

      @exception Exception::InvalidParameter is thrown if errors occur during the check
    */
    void checkDefaults(const String& name, const Param& defaults, const String& prefix = "") const;
    //@}

    ///@name Restriction handling
    //@{
    /**
      @brief Sets the valid strings for the parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::InvalidParameter is thrown, if one of the strings contains a comma character
      @exception Exception::ElementNotFound exception is thrown, if the parameter is no string parameter
    */
    void setValidStrings(const String& key, const std::vector<String>& strings);

    /**
      @brief Sets the minimum value for the integer or integer list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMinInt(const String& key, Int min);

    /**
      @brief Sets the maximum value for the integer or integer list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMaxInt(const String& key, Int max);

    /**
      @brief Sets the minimum value for the floating point or floating point list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMinFloat(const String& key, double min);

    /**
      @brief Sets the maximum value for the floating point or floating point list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMaxFloat(const String& key, double max);
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
      "prefix:misc" -> list("misc1","misc2")<BR>

      @param argc argc variable from command line
      @param argv argv variable from command line
      @param prefix prefix for all options
    */
    void parseCommandLine(const int argc, const char** argv, const String& prefix = "");

    /**
      @brief Parses command line arguments to specified key locations.

      Parses command line arguments to specified key locations and stores the result internally.

      @param argc argc variable from command line
      @param argv argv variable from command line
      @param options_with_one_argument a map of options that are followed by one argument (with key where they are stored)
      @param options_without_argument a map of options that are not followed by an argument (with key where they are stored). Options specified on the command line are set to the string 'true'.
      @param options_with_multiple_argument a map of options that are followed by several arguments (with key where they are stored)
      @param misc key where a StringList of all non-option arguments are stored
      @param unknown key where a StringList of all unknown options are stored
    */
    void parseCommandLine(const int argc, const char** argv, const Map<String, String>& options_with_one_argument, const Map<String, String>& options_without_argument, const Map<String, String>& options_with_multiple_argument, const String& misc = "misc", const String& unknown = "unknown");

    //@}

protected:

    /**
      @brief Returns a mutable reference to a parameter entry.

      @exception Exception::ElementNotFound is thrown for unset parameters
    */
    ParamEntry& getEntry_(const String& key) const;

    /// Constructor from a node which is used as root node
    Param(const Param::ParamNode& node);

    /// Invisible root node that stores all the data
    mutable Param::ParamNode root_;
  };

  /// Output of Param to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Param& param);

} // namespace OpenMS

