// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_

#include <seqan/arg_parse/arg_parse_exceptions.h>
#include <seqan/arg_parse/arg_parse_type_support.h>

#include <string>
#include <vector>

#include <sstream>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// ----------------------------------------------------------------------------
// Class ArgParseArgument
// ----------------------------------------------------------------------------

/*!
 * @class ArgParseArgument
 * @extends AssignableConcept
 * @headerfile <seqan/arg_parse.h>
 * @brief Information for a specific command line argument.
 *
 * @signature class ArgParseArgument
 */

/*!
 * @enum ArgParseArgument::ArgumentType
 * @headerfile <seqan/arg_parse.h>
 * @brief Define the type of an @link ArgParseArgument @endlink.
 *
 * @signature enum ArgParseArgument::ArgumentType;
 */

/*!
 * @var ArgParseArgument::ArgumentType ArgParseArgument::STRING
 * @brief Argument is a string.
 *
 * @var ArgParseArgument::ArgumentType ArgParseArgument::INTEGER
 * @brief Argument is a signed 32 bit integer.
 *
 * @var ArgParseArgument::ArgumentType ArgParseArgument::INT64
 * @brief Argument is a signed 64 bit integer.
 *
 * @var ArgParseArgument::ArgumentType ArgParseArgument::DOUBLE
 * @brief Argument is a floating point number stored as double.
 *
 * @var ArgParseArgument::ArgumentType ArgParseArgument::INPUTFILE
 * @brief Argument is an input file.
 *
 * @var ArgParseArgument::ArgumentType ArgParseArgument::OUTPUTFILE
 * @brief Argument is an output file.
 */

/**
.Class.ArgParseArgument
..cat:Miscellaneous
..summary:Stores information for a specific command line argument. It can be either an argument of
a ArgParseArgument or directly an Argument on the command line.
..signature:ArgParseArgument
..include:seqan/arg_parse.h
..see:Class.ArgParseOption
..see:Class.ArgumentParser
*/

/*!
 * @fn ArgParseArgument::ArgParseArgument
 * @brief Constructor
 *
 * @signature ArgParseArgument::ArgParseArgument(argumentType[, argumentLabel[, isListArgument[, numberOfArgument]]]);
 *
 * @param argumentType      Type of the argument (<tt>ArgParseArgument::ArgumentType</tt>).
 * @param argumentLabel     Label for the argument (<tt>char const *</tt>).
 * @param isListArgument    Whether or not this argument can be given multiple times (<tt>bool</tt>).
 * @param numberOfArguments Number of times the argument must be given.  E.g. set to 2 for the parser to always
 *                          expect two values (<tt>int</tt>, default is 1).
 */

/**
.Memfunc.ArgParseArgument#ArgParseArgument
..class:Class.ArgParseArgument
..summary:Constructor
..signature:ArgParseArgument (argumentType [, argumentLabel, isListArgument, numberOfArguments])
..param.argumentType:A ArgParseArgument.ArgumentType value defining the type (e.g., String) of the
ArgParseArgument.
...tableheader:Flag|Description
...table:$ArgParseArgument::STRING$|Argument is a string
...table:$ArgParseArgument::INTEGER$|Argument is an integer
...table:$ArgParseArgument::INT64|Argument is a 64 bit integer
...table:$ArgParseArgument::DOUBLE$|A float
...table:$ArgParseArgument::INPUTFILE$|An input file
...table:$ArgParseArgument::OUTPUTFILE$|An output file
 ..param.argumentLabel:Defines a user defined argument label for the help output. If this option is
 not set, ArgParseArgument will automatically define a label based on the ArgumentType.
..param.isListArgument:Defines if the argument can be given multiple times.
...default:false.
..param.numberOfArguments: Defines if the argument consists of defined number of elements (e.g., if
you want to provide an interval you would set this option to 2, so the parser knows that he needs
to search for exactly 2 values).
...default:1.
*/

class ArgParseArgument
{
public:
    enum ArgumentType
    {
        // argument is
        STRING,     // .. a string
        INTEGER,    // .. an integer
        INT64,      // .. a 64 bit integer
        DOUBLE,     // .. a float
        INPUTFILE,  // .. an inputfile (implicitly also a string)
        OUTPUTFILE  // .. an outputfile (implicitly also a string)
    };


    // ----------------------------------------------------------------------------
    // Members to store type information
    // ----------------------------------------------------------------------------
    ArgumentType _argumentType;
    unsigned     _numberOfValues;
    std::string  _argumentLabel;
    bool         _isListArgument;

    // ----------------------------------------------------------------------------
    // Members to store the values
    // ----------------------------------------------------------------------------
    std::vector<std::string>  defaultValue;
    std::vector<std::string>  value;

    // ----------------------------------------------------------------------------
    // Members for restrictions
    // ----------------------------------------------------------------------------
    std::string           minValue;
    std::string           maxValue;
    std::vector<std::string> validValues;

    // ----------------------------------------------------------------------------
    // Members to help text
    // ----------------------------------------------------------------------------
    std::string         _helpText;    // The help text shown on the command line

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------

    ArgParseArgument(ArgumentType argumentType,
                     std::string const & argumentLabel = "",
                     bool isListArgument = false,
                     unsigned numberOfValues = 1) :
        _argumentType(argumentType),
        _numberOfValues(numberOfValues),
        _argumentLabel(argumentLabel),
        _isListArgument(isListArgument),
        minValue(""),
        maxValue(""),
        _helpText("")
    {}

    ArgParseArgument(ArgParseArgument const & other) :
        _argumentType(other._argumentType),
        _numberOfValues(other._numberOfValues),
        _argumentLabel(other._argumentLabel),
        _isListArgument(other._isListArgument),
        defaultValue(other.defaultValue),
        value(other.value),
        minValue(other.minValue),
        maxValue(other.maxValue),
        validValues(other.validValues),
        _helpText(other._helpText)
    {}

    // TODO(holtgrew): This is superflous, generated by the compiler.
    ArgParseArgument & operator=(ArgParseArgument const & other)
    {
        if (this != &other)
        {
            _argumentType = other._argumentType;
            _numberOfValues = other._numberOfValues;
            _argumentLabel = other._argumentLabel;
            _isListArgument = other._isListArgument;
            defaultValue = other.defaultValue;
            value = other.value;
            minValue = other.minValue;
            maxValue = other.maxValue;
            validValues = other.validValues;
            _helpText = other._helpText;
        }

        return *this;
    }

};

// ==========================================================================
// Metafunctions
// ==========================================================================

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Helper Function _typeToString()
// ----------------------------------------------------------------------------

inline std::string _typeToString(ArgParseArgument const & me)
{
    std::string typeName = "";

    switch (me._argumentType)
    {
    case ArgParseArgument::DOUBLE:
        typeName = "double";
        break;

    case ArgParseArgument::INTEGER:
        typeName = "integer";
        break;

    case ArgParseArgument::INT64:
        typeName = "int64";
        break;

    case ArgParseArgument::STRING:
        typeName = "string";
        break;

    case ArgParseArgument::INPUTFILE:
        typeName = "inputfile";
        break;

    case ArgParseArgument::OUTPUTFILE:
        typeName = "outputfile";
        break;

    default:
        typeName = "unknown";
        break;
    }

    return typeName;
}

// ----------------------------------------------------------------------------
// Function isListArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isListArgument
 * @headerfile <seqan/arg_parse.h>
 *
 * @brief Returns whether the argument can be given more than one time.
 *
 * @signature bool isListArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it can be given multiple times, <tt>false</tt> otherwise.
 */

/**
.Function.isListArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument can be given multiple times.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument can be given multiple times.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.isListArgument
..include:seqan/arg_parse.h
*/

inline bool isListArgument(ArgParseArgument const & me)
{
    return me._isListArgument;
}

// ----------------------------------------------------------------------------
// Function isStringArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isStringArgument
 * @headerfile <seqan/arg_parse.h>
 *
 * @brief Returns whether the argument is a string.
 *
 * @signature bool ArgParseArgument#isStringArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a string, <tt>false</tt> otherwise.
 */

/**
.Function.isStringArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is a string.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument is a string argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isStringArgument(ArgParseArgument const & me)
{
    return (me._argumentType == ArgParseArgument::STRING) ||
           (me._argumentType == ArgParseArgument::INPUTFILE) ||
           (me._argumentType == ArgParseArgument::OUTPUTFILE);
}

// ----------------------------------------------------------------------------
// Function isIntegerArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isIntegerArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is an integer.
 *
 * @signature bool isIntegerArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is an integer, <tt>false</tt> otherwise.
 */

/**
.Function.isIntegerArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is an integer.
..cat:Miscellaneous
..signature:isIntegerArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument is an integer argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isIntegerArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INTEGER;
}

// ----------------------------------------------------------------------------
// Function isInt64Argument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isInt64Argument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a 64 bit integer.
 *
 * @signature bool isInt64Argument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a 64 bit integer, <tt>false</tt> otherwise.
 */

/**
.Function.isInt64Argument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is a 64 bit integer.
..cat:Miscellaneous
..signature:isInt64Argument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument is a 64 bit  integer argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isInt64Argument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INT64;
}

// ----------------------------------------------------------------------------
// Function isDoubleArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isDoubleArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a double integer.
 *
 * @signature bool isDoubleArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a double argument, <tt>false</tt> otherwise.
 */

/**
.Function.isDoubleArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is a double.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument is a double argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isDoubleArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::DOUBLE;
}

// ----------------------------------------------------------------------------
// Function isInputFileArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isInputFileArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a input file integer.
 *
 * @signature bool isInputFileArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a input file argument, <tt>false</tt> otherwise.
 */

/**
.Function.isInputFileArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is an input file.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if the argument argument is an input file argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isInputFileArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::INPUTFILE;
}

// ----------------------------------------------------------------------------
// Function isOutputFileArgument()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isOutputFileArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument is a output file integer.
 *
 * @signature bool isOutputFileArgument(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if it is a output file argument, <tt>false</tt> otherwise.
 */

/**
.Function.isOutputFileArgument
..class:Class.ArgParseArgument
..summary:Returns whether the argument is an output file.
..cat:Miscellaneous
..signature:isListArgument(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
...type:Class.ArgParseOption
..returns:$true$ if the argument argument is an output file argument.
..see:Memfunc.ArgParseArgument#ArgParseArgument.param.argumentType
..include:seqan/arg_parse.h
*/

inline bool isOutputFileArgument(ArgParseArgument const & me)
{
    return me._argumentType == ArgParseArgument::OUTPUTFILE;
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentLabel
 * @headerfile <seqan/arg_parse.h>
 * @brief Return argument label.
 *
 * @signature std::string getArgumentLabel(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return std::string The argument label as a STL string.
 */

/**
.Function.getArgumentLabel
..class:Class.ArgParseArgument
..summary:Returns the label for the given @Class.ArgParseArgument@. Either the user defined label
is returned or a default label (based on the ArgumentType is used).
..cat:Miscellaneous
..signature:getArgumentLabel(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:A $ShortCut.std::string$ containing the label.
..include:seqan/arg_parse.h
*/

inline std::string const getArgumentLabel(ArgParseArgument const & me)
{
    if (me._argumentLabel != "")
    {
        return me._argumentLabel;
    }
    else
    {
        // infer from argument type
        std::string baseLabel = "";
        if (isInputFileArgument(me) || isOutputFileArgument(me))
            baseLabel = "FILE";
        else if (isStringArgument(me))
            baseLabel = "STR";
        else if (isIntegerArgument(me) || isDoubleArgument(me))
            baseLabel = "NUM";

        std::string finalLabel;

        if (me._numberOfValues != 1)
        {
            for (unsigned i = 0; i < me._numberOfValues; ++i)
            {
                if (i != 0)
                    append(finalLabel, " ");
                append(finalLabel, baseLabel);
            }
        }
        else if (isListArgument(me))
            finalLabel = baseLabel;                         // maybe we want to customize list labels
        else
            finalLabel = baseLabel;

        return finalLabel;
    }
}

// ----------------------------------------------------------------------------
// Helper Function _intervalAssert()
// ----------------------------------------------------------------------------

// this methods ensures that the given arguments define a non emtpy value interval
// otherwise it will trigger a SEQAN_CHECK failure
template <typename TIntervalBorder>
inline void _intervalAssert(const std::string minValueAsString, const std::string maxValueAsString)
{
    if (minValueAsString != "" && maxValueAsString != "")
        SEQAN_CHECK(_cast<TIntervalBorder>(minValueAsString) < _cast<TIntervalBorder>(maxValueAsString),
                    "The interval [%s:%s] is empty. Please specify a valid, non-empty interval.",
                    minValueAsString.c_str(),
                    maxValueAsString.c_str());
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setMinValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set smallest allowed value for the argument.
 *
 * @signature void setMinValue(arg, minValue);
 *
 * @param arg      The ArgParseArgument to set the smallest value of.
 * @param minValue The smallest value to set (<tt>std::string</tt>).
 *
 * @return std::string The argument label as a STL string.
 */

/**
.Function.setMinValue
..class:Class.ArgParseArgument
..summary:Sets the minimum value of a @Class.ArgParseArgument@ object.
..cat:Miscellaneous
..signature:setMinValue(argument,minValue)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.minValue:A std::string containing a string representation of the minimum value
of the @Class.ArgParseArgument@.
..include:seqan/arg_parse.h
*/

inline void setMinValue(ArgParseArgument & me, const std::string minValue)
{
    if (isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(minValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else if (isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(minValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else if (isInt64Argument(me))
    {
        SEQAN_CHECK(_isCastable<__int64>(minValue), "The maximal value for a 64 integer argument must be a 64 bit integer");
        _intervalAssert<__int64>(minValue, me.maxValue);
        me.minValue = minValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setMaxValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set smallest allowed value for the argument.
 *
 * @signature void setMaxValue(arg, maxValue);
 *
 * @param arg      The ArgParseArgument to set the smallest value of.
 * @param maxValue The largest value to set (<tt>std::string</tt>).
 *
 * @return std::string The argument label as a STL string.
 */

/**
.Function.setMaxValue
..class:Class.ArgParseArgument
..summary:Sets the maximum value of a @Class.ArgParseArgument@ object.
..cat:Miscellaneous
..signature:setMaxValue(argument,maxValue)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.maxValue:A std::string containing a string representation of the maximum value
of the @Class.ArgParseArgument@.
..include:seqan/arg_parse.h
*/

inline void setMaxValue(ArgParseArgument & me, const std::string maxValue)
{
    if (isDoubleArgument(me))
    {
        SEQAN_CHECK(_isCastable<double>(maxValue), "The maximal value for a double argument must be double.");
        _intervalAssert<double>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else if (isIntegerArgument(me))
    {
        SEQAN_CHECK(_isCastable<int>(maxValue), "The maximal value for an integer argument must be an integer");
        _intervalAssert<int>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else if (isInt64Argument(me))
    {
        SEQAN_CHECK(_isCastable<int>(maxValue), "The maximal value for a 64 bit integer argument must be an 64 bit integer");
        _intervalAssert<int>(me.minValue, maxValue);
        me.maxValue = maxValue;
    }
    else
        SEQAN_FAIL("min/max values are not applicable to non numeric arguments");
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setValidValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Set list of valid values.
 *
 * @signature void setValidValues(arg, values);
 *
 * @param arg    The ArgParseArgument to set the valid values for.
 * @param values Either a <tt>std::string</tt> containing all valid entries, separated by spaces or a
 *               <tt>std::vector&lt;std::string&gt;</tt> with the valid entries.
 *
 * If the argument is of type string then the list of valid values is the case-sensitive list of string values
 * allowed for this argument.  If it is an input or output file then the list of valid values is a list of
 * case-insentive file extensions identifying the allowed types.
 *
 * @section Examples
 *
 * An example of setting allowed values for a string option.
 *
 * @code{.cpp}
 * seqan::ArgParseArgument stringArg(seqan::ArgParseArgument::STRING);
 * setValidValues(stringArg, "one two three");  // one of {"one", "two", "three"}
 *
 * std::vector<std::string> values;
 * values.push_back("four");
 * values.push_back("five");
 * setValidValues(stringArg, values);  // one of {"four", "five"}
 * @endcode
 *
 * An example for an input file option.  Note that by changing <tt>INPUTFILE</tt> to <tt>OUTPUTFILE</tt> below,
 * the example would be the same for output files.
 *
 * @code{.cpp}
 * seqan::ArgParseArgument fileArg(seqan::ArgParseArgument::INPUTFILE);
 * setValidValues(fileArg, "fq fastq");  // file must end in ".fq" or ".fastq"
 *
 * std::vector<std::string> values;
 * values.push_back("sam");
 * values.push_back("bam");
 * setValidValues(fileArg, values);  // file must end in ".sam" or ".bam"
 * @endcode
 */

/**
.Function.setValidValues
..class:Class.ArgParseArgument
..summary:Sets the set of allowed values of a @Class.ArgParseArgument@ object.
..cat:Miscellaneous
..signature:setValidValues(argument,values)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.values:A std::vector<std::string> containing all valid entries for the option or a std::string
with valid values separated by spaces.
..remarks:If the argument or option is an in- or output file. The valid strings will be interpreted as
file endings and the command line parser checks if the provided file has the required file ending.
..include:seqan/arg_parse.h
*/

inline void setValidValues(ArgParseArgument & me, std::vector<std::string> const & values)
{
    if (isDoubleArgument(me) || isIntegerArgument(me))
        SEQAN_FAIL("ArgParseArgument does not support setting valid values for numeric arguments.");

    me.validValues = values;
}

inline void setValidValues(ArgParseArgument & me, std::string const & valuesString)
{
    // convert array to String<std::string>
    std::vector<std::string> values;
    std::string current_argument;

    for (std::string::const_iterator ch  = valuesString.begin(); ch != valuesString.end(); ++ch)
    {
        if (*ch == ' ')
        {
            appendValue(values, current_argument);
            current_argument = "";
        }
        else
        {
            append(current_argument, *ch);
        }
    }
    if (current_argument != "")
        appendValue(values, current_argument);

    setValidValues(me, values);
}

// ----------------------------------------------------------------------------
// Function setHelpText()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#setHelpText
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the help text for an ArgParseArgument.
 *
 * @signature void setHelpText(arg, text);
 *
 * @param arg  The ArgParseArgument to set the help text for.
 * @param text The text to display as the description of the argument (<tt>std::string</tt>).
 */

/**
.Function.setHelpText
..class:Class.ArgParseArgument
..summary:Sets the help text for an ArgParseArgument.
..cat:Miscellaneous
..signature:setHelpText(argument,text)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.text:A std::string describing the argument.
..include:seqan/arg_parse.h
*/

inline void setHelpText(ArgParseArgument & me, std::string const & text)
{
    me._helpText = text;
}

// ----------------------------------------------------------------------------
// Helper Function _isInInterval()
// ----------------------------------------------------------------------------

// check if the given value is in the provided interval
template <typename TTarget, typename TString>
inline bool _isInInterval(TString value, TString lowerIntervalBound, TString upperIntervalBound)
{
    bool isInInterval = true;

    if (lowerIntervalBound != "")
        isInInterval &= (_cast<TTarget>(lowerIntervalBound) <= _cast<TTarget>(value));
    if (upperIntervalBound != "")
        isInInterval &= (_cast<TTarget>(value) <= _cast<TTarget>(upperIntervalBound));

    return isInInterval;
}

// ----------------------------------------------------------------------------
// Helper Function _checkNumericArgument()
// ----------------------------------------------------------------------------

// test if the values can be assigned to the option and is in the given boundaries
template <typename TNumerical>
inline void _checkNumericArgument(ArgParseArgument const & me, std::string const & value)
{
    if (!_isCastable<TNumerical>(value))
    {
        std::stringstream what;
        what << "the given value '" << value << "' cannot be casted to " << _typeToString(me);
        throw ParseException(what.str());
    }

    if (!_isInInterval<TNumerical>(value, me.minValue, me.maxValue))
    {
        std::stringstream what;
        what << "the given value '" << value << "' is not in the interval ["
             << (me.minValue != "" ? me.minValue : "-inf") << ":"
             << (me.maxValue != "" ? me.maxValue : "+inf") << "]";

        throw ParseException(what.str());
    }
}

// ----------------------------------------------------------------------------
// Helper Function _compareExtension()
// ----------------------------------------------------------------------------

inline bool _compareExtension(std::string const & str, std::string const & ext)
{
    std::string str_ext = str.substr(str.size() - ext.size());
    for (size_t i = 0; i < str_ext.size() && i < ext.size(); ++i)
    {
        if (tolower(str_ext[i]) != tolower(ext[i]))
            return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Helper Function _checkStringRestrictions()
// ----------------------------------------------------------------------------

inline void _checkStringRestrictions(ArgParseArgument const & me, std::string const & value)
{
    typedef std::vector<std::string>::const_iterator TVectorIterator;

    if (!empty(me.validValues))
    {
        bool isContained = false;
        for (TVectorIterator validValue = me.validValues.begin();
             validValue != me.validValues.end();
             ++validValue)
        {
            // if it is an input or output file, we only check the file endings
            if (isInputFileArgument(me) || isOutputFileArgument(me))
            {
                if (length(*validValue) > length(value))
                    continue;
                else
                    isContained |= _compareExtension(value, *validValue);
            }
            else
            {
                isContained |= (*validValue == value);
            }
            if (isContained)
                break;
        }
        if (!isContained)
        {
            std::stringstream what;
            what << "the given value '" << value <<
            "' is not in the list of allowed" <<
            ((isInputFileArgument(me) || isOutputFileArgument(me)) ? " file extensions " : " values ") <<
            "[";
            for (TVectorIterator validValue = me.validValues.begin();
                 validValue != me.validValues.end();
                 ++validValue)
            {
                if (validValue != me.validValues.begin())
                    what << ", ";
                what << ((isInputFileArgument(me) || isOutputFileArgument(me)) ? "*" : "") << *validValue;
            }
            what << "]";
            throw ParseException(what.str());
        }
    }
}

// ----------------------------------------------------------------------------
// Function _checkValue()
// ----------------------------------------------------------------------------

inline void _checkValue(ArgParseArgument const & me, std::string const & value)
{
    // type checks
    if (isIntegerArgument(me))
        _checkNumericArgument<int>(me, value);

    if (isInt64Argument(me))
        _checkNumericArgument<__int64>(me, value);

    if (isDoubleArgument(me))
        _checkNumericArgument<double>(me, value);

    // check valid values
    if (isStringArgument(me))
        _checkStringRestrictions(me, value);
}

// ----------------------------------------------------------------------------
// Function _assignArgumentValue()
// ----------------------------------------------------------------------------

/**
.Internal.Function._assignArgumentValue
..class:Class.ArgParseArgument
..summary:Assigns the given value (if applicable) to the @Class.ArgParseArgument@ object. If
the @Class.ArgParseArgument@ is a list or can hold multiple values
(@Memfunc.ArgParseArgument#ArgParseArgument.param.numberOfArguments@) the value will be appended.
Otherwise the value will be overwritten.
..cat:internal
..signature:_assignArgumentValue(argument,value [, argNo])
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.value:A std::string containing the value that should be assigned.
..include:seqan/arg_parse.h
*/

inline void _assignArgumentValue(ArgParseArgument & me, std::string const & value) throw (ParseException)
{
    // check values
    _checkValue(me, value);

    // assignment
    if (isListArgument(me)) // just append
        appendValue(me.value, value, Exact());
    else
    {
        // check if we already set all expected arguments
        if (length(me.value) == me._numberOfValues)
            clear(me.value);
        appendValue(me.value, value, Exact());
    }
}

// ----------------------------------------------------------------------------
// Function getArgumentValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Return the value of the argument.
 *
 * @signature std::string getArgumentValue(arg[, argNo]);
 *
 * @param arg   The ArgParseArgument to query.
 * @param argNo In case that the ArgParseArgument allowed multiple values, give the index of the argument
 *              that you want to retrieve (<tt>unsigned</tt>, starts at 0).
 *
 * @return std::string Const-reference to the argument value.
 */

/**
.Function.ArgParseArgument#getArgumentValue
..class:Class.ArgParseArgument
..summary:Returns the value of the @Class.ArgParseArgument@ object. If
the @Class.ArgParseArgument@ is a list or can hold multiple values
(@Memfunc.ArgParseArgument#ArgParseArgument.param.numberOfArguments@) you can specify which value
you want to get. If not set the first value will be returned.
..cat:Miscellaneous
..signature:getArgumentValue(argument [, argNo])
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.argNo:If the argument is a list,  the $argNo$-th list element is returned.
..returns:The value set at position $position$.
..include:seqan/arg_parse.h
*/

inline std::string const & getArgumentValue(ArgParseArgument const & me, unsigned argNo)
{
    SEQAN_CHECK(argNo < me.value.size() || argNo < me.defaultValue.size(),
                "ArgParseArgument: No value set for index %d", argNo);

    if (argNo < me.value.size())
        return me.value[argNo];
    else
        return me.defaultValue[argNo];
}

inline std::string const & getArgumentValue(ArgParseArgument const & me)
{
    return getArgumentValue(me, 0);
}

// ----------------------------------------------------------------------------
// Function getArgumentValues()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#getArgumentValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Return all values of the argument.
 *
 * @signature std::vector<std::string> getArgumentValue(arg);
 *
 * @param arg   The ArgParseArgument to query.
 *
 * @return std::vector<std::string> Const-reference to the argument values.
 */

/**
.Function.getArgumentValues
..class:Class.ArgParseArgument
..summary:Returns all values of the @Class.ArgParseArgument@ object as const std::vector.
..cat:Miscellaneous
..signature:getArgumentValues(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$std::vector<std::string>$ containing the values. If no value was set and no
default value exists an empty vector will be returned.
..include:seqan/arg_parse.h
*/

inline std::vector<std::string> const & getArgumentValues(ArgParseArgument const & me)
{
    if (!me.value.empty())
        return me.value;
    else
        return me.defaultValue;
}

// ----------------------------------------------------------------------------
// Function hasArgumentValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#hasArgumentValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Return whether a value is available.
 *
 * @signature bool hasValue(arg[, pos]);
 *
 * @param arg The ArgParseArgument to query.
 * @param pos The position of the argument in case of being a list (<tt>unsigned</tt>, 0-based, default is 0).
 *
 * @return bool <tt>true</tt> if <tt>pos</tt> is less than the size and the argument is non-empty.
 */

/**
.Function.ArgParseArgument#hasValue
..class:Class.ArgParseArgument
..summary:Returns true if a value for the given position is available.
..cat:Miscellaneous
..signature:hasValue(argument [, position=0])
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..param.position:The position for which the availability should be tested.
..returns: $true$ if a value is available, $false$ if not.
..include:seqan/arg_parse.h
*/
inline bool hasValue(ArgParseArgument const & arg, unsigned position)
{
    return arg.value.size() > position || arg.defaultValue.size() > position;
}

inline bool hasValue(ArgParseArgument const & arg)
{
    return hasValue(arg, 0);
}

// ----------------------------------------------------------------------------
// Function isSet()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#isSet
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns true if a value was assigned to the argument.
 *
 * @signature bool isSet(arg):
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return bool <tt>true</tt> if a value was assigned, <tt>false</tt> otherwise.
 */

/**
.Function.ArgParseArgument#isSet
..class:Class.ArgParseArgument
..summary:Returns true if a value was assigned to the argument.
..cat:Miscellaneous
..signature:isSet(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if a value was assigned to the argument, $false$ if not.
..include:seqan/arg_parse.h
*/

inline bool isSet(ArgParseArgument const & me)
{
    return !me.value.empty();
}

// ----------------------------------------------------------------------------
// Function hasDefault()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#hasDefault
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the argument has a default value.
 *
 * @signature bool hasDefault(arg);
 *
 * @param arg The argument to query.
 *
 * @return bool <tt>true</tt> if the argument has a default value and <tt>false</tt> if not.
 */

/**
.Function.ArgParseArgument#hasDefault
..summary:Returns true if a default value was given for that argument.
..cat:Miscellaneous
..signature:hasDefault(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:$true$ if a default value was given for the argument, $false$ if not.
..include:seqan/arg_parse.h
*/

inline bool hasDefault(ArgParseArgument const & me)
{
    return !me.defaultValue.empty();
}

// ----------------------------------------------------------------------------
// Function numberOfArguments
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseArgument#numberOfAllowedValues
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns the number of allowed values for this ArgParseArgument.
 *
 * @signature unsigned numberOfAllowedValues(arg);
 *
 * @param arg The ArgParseArgument to query.
 *
 * @return unsigned The number of allowed values.
 */

/**
.Function.numberOfAllowedValues
..class:Class.ArgParseArgument
..summary:Returns the number of allowed values for this @Class.ArgParseArgument@.
..cat:Miscellaneous
..signature:numberOfAllowedValues(argument)
..param.argument:The @Class.ArgParseArgument@ object.
...type:Class.ArgParseArgument
..returns:The number of allowed values for this @Class.ArgParseArgument@.
..include:seqan/arg_parse.h
*/

inline unsigned numberOfAllowedValues(ArgParseArgument const & me)
{
    return me._numberOfValues;
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_ARGUMENT_H_
