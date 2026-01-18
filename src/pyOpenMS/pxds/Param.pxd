from libcpp.string cimport string as libcpp_utf8_string
from libcpp.string cimport string as libcpp_string
from Types cimport *
from ParamValue cimport *
from ParamNode cimport *
from ParamEntry cimport *


# this class has addons, see the ./addons folder

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS":

    # pythonic helper functions in ../addons/Param.pyx !!!!

    cdef cppclass Param:
         # wrap-doc:
         #  Management and storage of parameters / INI files.
         #  This class provides a means to associate string names to int/double/string/StringList values.
         #  It allows for parameter hierarchies and to save/load the data as XML.
         #  Hierarchy levels are separated from each other by colons (e.g., 'algorithm:common:threshold').
         #  Each parameter and section has a description. Newline characters in the description are possible.
         #  Each parameter can be annotated with an arbitrary number of tags (e.g., 'advanced').

         # COMMENT: Helper functions for Python (Pythonic dict-like interface)
         from_dict(d) # wrap-ignore
         to_dict() # wrap-ignore
         initPluginParam(name, version) # wrap-ignore
         asDict() # wrap-ignore
         keys() # wrap-ignore
         items() # wrap-ignore
         values() # wrap-ignore
         update(dict) # wrap-ignore
         get(key, default=None) # wrap-ignore
         __getitem__(key) # wrap-ignore
         __setitem__(key, value) # wrap-ignore
         __iter__() # wrap-ignore
         __len__() # wrap-ignore
         __contains__(key) # wrap-ignore

         Param() except + nogil  # wrap-doc:Default constructor
         Param(Param &) except + nogil  # wrap-doc:Copy constructor
         bool operator==(Param) except + nogil  # wrap-doc:Equality operator

         # wrap-doc:
         #  Sets a value. The key can contain colons to specify sections (e.g., 'section:subsection:parameter').
         #  Optionally, a description and tags can be provided
         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc, libcpp_vector[libcpp_utf8_string] tags) except + nogil
         void setValue(libcpp_utf8_string key, ParamValue val, libcpp_utf8_string desc) except + nogil
         void setValue(libcpp_utf8_string key, ParamValue val) except + nogil
         ParamValue getValue(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns the value of the parameter specified by key. Raises exception if not found
         ValueType getValueType(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns the type of the parameter specified by key. Raises exception if not found

         ParamEntry getEntry(libcpp_utf8_string) except + nogil  # wrap-doc:Returns the whole parameter entry (value, description, tags, restrictions). Raises exception if not found
         bool exists(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns True if the parameter exists, False otherwise

         void addTag(libcpp_utf8_string key, libcpp_utf8_string tag) except + nogil  # wrap-doc:Adds a tag to the entry specified by key (e.g., 'advanced', 'required', 'input file')
         void addTags(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] tags) except + nogil  # wrap-doc:Adds multiple tags to the entry specified by key
         int hasTag(libcpp_utf8_string key, libcpp_utf8_string tag) except + nogil  # wrap-doc:Returns True if the parameter has the specified tag
         libcpp_vector[libcpp_string] getTags(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns all tags of the entry specified by key
         void clearTags(libcpp_utf8_string key) except + nogil  # wrap-doc:Removes all tags from the entry specified by key

         libcpp_utf8_output_string getDescription(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns the description of the parameter specified by key
         void setSectionDescription(libcpp_utf8_string key, libcpp_utf8_string desc) except + nogil  # wrap-doc:Sets a description for an existing section (not for values)
         libcpp_utf8_output_string getSectionDescription(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns the description of the section specified by key (empty string if not found)

         void addSection(libcpp_utf8_string key, libcpp_utf8_string desc) except + nogil  # wrap-doc:Adds a parameter section with the given key and description

         Size size() except + nogil  # wrap-doc:Returns the number of parameter entries (leaves)
         bool empty() except + nogil  # wrap-doc:Returns True if there are no entries

         void clear() except + nogil  # wrap-doc:Deletes all entries
         void insert(libcpp_utf8_string prefix, Param param) except + nogil  # wrap-doc:Inserts all values of another Param object with the given prefix

         void remove(libcpp_utf8_string key) except + nogil  # wrap-doc:Removes an entry or section (when key ends with ':') by exact name match
         void removeAll(libcpp_utf8_string prefix) except + nogil  # wrap-doc:Removes all entries and sections that start with the given prefix

         # wrap-doc:
         #  Returns a new Param containing all entries that start with the given prefix.
         #  If remove_prefix is True, the prefix is removed from the keys in the returned Param
         Param copy(libcpp_utf8_string prefix, bool) except + nogil
         Param copy(libcpp_utf8_string prefix) except + nogil

         # wrapped manually for overloading with dict parameter:
         bool update(Param p_old, bool add_unknow) except + nogil  # wrap-ignore
         bool update(Param p_old) except + nogil  # wrap-ignore

         void merge(Param toMerge) except + nogil  # wrap-doc:Adds missing parameters from another Param object without modifying existing ones

         # wrap-doc:
         #  Inserts all values from defaults that are not already set.
         #  Optionally adds a prefix to all keys and prints a message for each default value set
         void setDefaults(Param defaults, libcpp_utf8_string prefix, bool showMessage) except + nogil
         void setDefaults(Param defaults, libcpp_utf8_string prefix) except + nogil
         void setDefaults(Param defaults) except + nogil

         # wrap-doc:
         #  Checks current parameter entries against given defaults.
         #  Validates types, string restrictions, and numeric ranges. Raises exception on invalid parameters
         void checkDefaults(libcpp_utf8_string name, Param defaults, libcpp_utf8_string prefix) except + nogil
         void checkDefaults(libcpp_utf8_string name, Param defaults) except + nogil

         libcpp_vector[libcpp_utf8_string] getValidStrings(libcpp_utf8_string key) except + nogil  # wrap-doc:Returns the list of valid string values for the parameter

         void setValidStrings(libcpp_utf8_string key, libcpp_vector[libcpp_utf8_string] strings) except + nogil  # wrap-doc:Sets the list of valid string values for the parameter (checked by checkDefaults)
         void setMinInt(libcpp_utf8_string key, int min) except + nogil  # wrap-doc:Sets the minimum allowed value for an integer parameter
         void setMaxInt(libcpp_utf8_string key, int max) except + nogil  # wrap-doc:Sets the maximum allowed value for an integer parameter
         void setMinFloat(libcpp_utf8_string key, double min) except + nogil  # wrap-doc:Sets the minimum allowed value for a float parameter
         void setMaxFloat(libcpp_utf8_string key, double max) except + nogil  # wrap-doc:Sets the maximum allowed value for a float parameter

         #void parseCommandLine(int argc, char ** argv, String prefix) # wrap-ignore
         #void parseCommandLine(int argc, char ** argv) # wrap-ignore

         ParamIterator begin() except + nogil  # wrap-ignore
         ParamIterator end()   except + nogil  # wrap-ignore

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param":

    cdef cppclass ParamIterator:
        # wrap-ignore
        # no-pxd-import
        ParamIterator operator++() except + nogil
        ParamIterator operator--() except + nogil
        libcpp_utf8_string getName() except + nogil
        int operator==(ParamIterator) except + nogil
        int operator!=(ParamIterator) except + nogil
        int operator<(ParamIterator) except + nogil
        int operator>(ParamIterator) except + nogil
        int operator<=(ParamIterator) except + nogil
        int operator>=(ParamIterator) except + nogil

        # Returns the traceback of the opened and closed sections
        libcpp_vector[TraceInfo] getTrace() except + nogil

cdef extern from "<OpenMS/DATASTRUCTURES/Param.h>" namespace "OpenMS::Param::ParamIterator":

    cdef cppclass TraceInfo:
        # wrap-doc:
        #  Struct that captures information on entered/left nodes during Param iteration.
        #  Used to track the traceback of opened and closed sections.

        TraceInfo(libcpp_utf8_string n, libcpp_utf8_string d, bool o) except + nogil  # wrap-doc:Constructor with name, description, and open flag
        TraceInfo(TraceInfo) except + nogil  # wrap-doc:Copy constructor

        # name of the node
        libcpp_string name
        # description of the node
        libcpp_string description
        # If it was opened (true) or closed (false)
        bool opened
