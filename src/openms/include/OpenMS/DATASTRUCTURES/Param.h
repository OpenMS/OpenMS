// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/ParamValue.h>
#include <OpenMS/OpenMSConfig.h>

#include <set>
#include <string>
#include <map>

namespace OpenMS
{

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
      ParamEntry(const std::string& n, const ParamValue& v, const std::string& d, const std::vector<std::string>& t = std::vector<std::string>());
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
      bool isValid(std::string& message) const;
      /// Equality operator (only name and value are compared)
      bool operator==(const ParamEntry& rhs) const;

      /// Name of the entry
      std::string name;
      /// Description of the entry
      std::string description;
      /// Value associated with the entry
      ParamValue value;
      /// Tags list, used e.g. for advanced parameter tag
      std::set<std::string> tags;
      ///@name Restrictions to accepted values (used in checkDefaults)
      //@{
      double min_float; ///< Default: - std::numeric_limits<double>::max()
      double max_float; ///< Default: std::numeric_limits<double>::max()
      int min_int; ///< Default: - std::numeric_limits<Int>::max()
      int max_int; ///< Default: std::numeric_limits<Int>::max()
      std::vector<std::string> valid_strings; ///< Default: empty
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
      ParamNode(const std::string& n, const std::string& d);
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
      EntryIterator findEntry(const std::string& name);
      /**
        @brief Look up subnode of this node (local search)

        Returns the end iterator if no entry is found
      */
      NodeIterator findNode(const std::string& name);
      /**
        @brief Look up the parent node of the entry or node corresponding to @p name (tree search)

        Returns 0 if no entry is found
      */
      ParamNode* findParentOf(const std::string& name);
      /**
        @brief Look up the entry corresponding to @p name (tree search)

        Returns 0 if no entry is found
      */
      ParamEntry* findEntryRecursive(const std::string& name);

      ///Inserts a @p node with the given @p prefix
      void insert(const ParamNode& node, const std::string& prefix = "");
      ///Inserts an @p entry with the given @p prefix
      void insert(const ParamEntry& entry, const std::string& prefix = "");
      ///Returns the number of entries in the whole subtree
      size_t size() const;
      ///Returns the name suffix of a @p key (the part behind the last ':' character)
      std::string suffix(const std::string& key) const;

      /// Name of the node
      std::string name;
      /// Description of the node
      std::string description;
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
        inline TraceInfo(const std::string& n, const std::string& d, bool o) :
          name(n),
          description(d),
          opened(o)
        {
        }

        /// name of the node
        std::string name;
        /// description of the node
        std::string description;
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
      std::string getName() const;
      /// Returns the traceback of the opened and closed sections
      const std::vector<TraceInfo>& getTrace() const;

protected:
      /// Pointer to the root node
      const Param::ParamNode* root_;
      /// Index of the current ParamEntry (-1 means invalid)
      int current_;
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
    void setValue(const std::string& key, const ParamValue& value, const std::string& description = "", const std::vector<std::string>& tags = std::vector<std::string>());

    /**
      @brief Returns a value of a parameter.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    const ParamValue& getValue(const std::string& key) const;

    /**
      @brief Returns the type of a parameter.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    ParamValue::ValueType getValueType(const std::string& key) const;

    /**
      @brief Returns the whole parameter entry.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    const ParamEntry& getEntry(const std::string& key) const;

    /**
      @brief Tests if a parameter is set (expecting its fully qualified name, e.g., TextExporter:1:proteins_only)

      @param key The fully qualified name of the parameter to check.
      @return True if the parameter exists, false otherwise.
    */
    bool exists(const std::string& key) const;

    /**
      @brief Checks whether a section is present.

      @param key The key of the section to be searched for. May or may not contain ":" suffix.
      @return True if the section exists, false otherwise.
     */
    bool hasSection(const std::string& key) const;

    /**
      @brief Find leaf node by name (if it exists).

      @param leaf The name of the parameter to find excluding the path parameter, e.g., given the parameter TextExporter:1:proteins_only the leaf would be named proteins_only.
      @return Returns end() if leaf does not exist.
    */
    ParamIterator findFirst(const std::string& leaf) const;

    /**
      @brief Find next leaf node by name (if it exists), not considering the @p start_leaf

      @param leaf The name of the parameter to find excluding the path parameter, e.g., given the parameter TextExporter:1:proteins_only the leaf would be named proteins_only.
      @param start_leaf The already found leaf, that should not be considered during this search.
      @return Returns end() if leaf does not exist.
    */
    ParamIterator findNext(const std::string& leaf, const ParamIterator& start_leaf) const;
    //@}

    ///@name Tags handling
    //@{

    /**
      @brief Adds the tag @p tag to the entry @p key

      E.g. "advanced", "required", "input file", "output file"

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
      @exception Exception::InvalidValue is thrown if the tag contain a comma character.
    */
    void addTag(const std::string& key, const std::string& tag);

    /**
      @brief Adds the tags in the list @p tags to the entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
      @exception Exception::InvalidValue is thrown if a tag contain a comma character.
    */
    void addTags(const std::string& key, const std::vector<std::string>& tags);

    /**
      @brief Returns if the parameter @p key has a tag

      Example: The tag 'advanced' is used in the GUI to determine which parameters are always displayed
      and which parameters are displayed only in 'advanced mode'.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    bool hasTag(const std::string& key, const std::string& tag) const;

    /**
      @brief Returns the tags of entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    std::vector<std::string> getTags(const std::string& key) const;

    /**
      @brief Removes all tags from the entry @p key

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    void clearTags(const std::string& key);
    //@}


    ///@name Descriptions handling
    //@{

    /**
      @brief Returns the description of a parameter.

      @exception Exception::ElementNotFound is thrown if the parameter does not exists.
    */
    const std::string& getDescription(const std::string& key) const;

    /**
      @brief Sets a description for an existing section

      Descriptions for values cannot be set with this method.
      They have to be set when inserting the value itself.

      @exception Exception::ElementNotFound is thrown if the section does not exists.
    */
    void setSectionDescription(const std::string& key, const std::string& description);

    /**
      @brief Returns the description corresponding to the section with name @p key.

      If the section does not exist an empty string is returned.
    */
    const std::string& getSectionDescription(const std::string& key) const;

    /**
    @brief Adds a parameter section under the path @p key with the given @p description.

    If the section already exists, the description is only overwritten if not empty.
    */
    void addSection(const std::string& key, const std::string& description);
    //@}

    ///@name Manipulation of the whole parameter set
    //@{

    ///Returns the number of entries (leafs).
    size_t size() const;

    ///Returns if there are no entries.
    bool empty() const;

    /// Deletes all entries
    void clear();

    /// Insert all values of @p param and adds the prefix @p prefix.
    /// You should append ':' to prefix manually when you want it to be a section.
    void insert(const std::string& prefix, const Param& param);

    /**
      @brief Remove the entry @p key or a section @p key (when suffix is ':')

      Remove deletes either an entry or a section (when @p key ends with ':'),
      by matching the exact name. No partial matches are accepted.

      If an empty internal node remains, the tree is pruned until every node has either a successor node
      or a leaf, i.e. no naked nodes remain.
    */
    void remove(const std::string& key);

    /**
      @brief Remove all entries that start with @p prefix

      Partial are valid as well. All entries and sections which match the prefix are deleted.

      If an empty internal node remains, the tree is pruned until every node has either a successor node
      or a leaf, i.e. no naked nodes remain.
    */
    void removeAll(const std::string& prefix);

    /**
      @brief Returns a new Param object containing all entries that start with @p prefix.

      @param prefix should contain a ':' at the end if you want to extract a subtree.
             Otherwise not only nodes, but as well values with that prefix are copied.
      @param remove_prefix indicates if the prefix is removed before adding entries to the new Param
    */
    Param copy(const std::string& prefix, bool remove_prefix = false) const;

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

      Calls update(p_outdated, true, add_unknown, false, false, OPENMS_LOG_WARN) and returns its value.
    */
    bool update(const Param& p_outdated, const bool add_unknown = false);

    /**
      @brief Rescue parameter <b>values</b> from @p p_outdated to current param

      Calls update(p_outdated, true, add_unknown, false, false, stream) and returns its value.
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
      @param add_unknown Add unknown parameters from @p p_outdated to this param object.
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
    void setDefaults(const Param& defaults, const std::string& prefix = "", bool showMessage = false);

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
    void checkDefaults(const std::string& name, const Param& defaults, const std::string& prefix = "") const;
    //@}

    ///@name Restriction handling
    //@{
    /**
      @brief Sets the valid strings for the parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::InvalidParameter is thrown, if one of the strings contains a comma character
      @exception Exception::ElementNotFound exception is thrown, if the parameter is no string/stringlist parameter
    */
    void setValidStrings(const std::string& key, const std::vector<std::string>& strings);

    /**
      @brief Gets he valid strings for the parameter @p key.

      @exception Exception::ElementNotFound exception is thrown, if the parameter is no string/stringlist parameter
    */
    const std::vector<std::string>& getValidStrings(const std::string& key) const;

    /**
      @brief Sets the minimum value for the integer or integer list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMinInt(const std::string& key, int min);

    /**
      @brief Sets the maximum value for the integer or integer list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMaxInt(const std::string& key, int max);

    /**
      @brief Sets the minimum value for the floating point or floating point list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMinFloat(const std::string& key, double min);

    /**
      @brief Sets the maximum value for the floating point or floating point list parameter @p key.

      It is only checked in checkDefaults().

      @exception Exception::ElementNotFound is thrown if @p key is not found or if the parameter type is wrong
    */
    void setMaxFloat(const std::string& key, double max);
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
    void parseCommandLine(const int argc, const char** argv, const std::string& prefix = "");

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
    void parseCommandLine(const int argc, const char** argv, const std::map<std::string, std::string>& options_with_one_argument, const std::map<std::string, std::string>& options_without_argument, const std::map<std::string, std::string>& options_with_multiple_argument, const std::string& misc = "misc", const std::string& unknown = "unknown");

    //@}

protected:

    /**
      @brief Returns a mutable reference to a parameter entry.

      @exception Exception::ElementNotFound is thrown for unset parameters
    */
    ParamEntry& getEntry_(const std::string& key) const;

    /// Constructor from a node which is used as root node
    Param(const Param::ParamNode& node);

    /// Invisible root node that stores all the data
    mutable Param::ParamNode root_;
  };

  /// Output of Param to a stream.
  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const Param& param);

} // namespace OpenMS

