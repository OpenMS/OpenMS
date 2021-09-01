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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_

#include <string>
#include <vector>
#include <seqan/arg_parse/arg_parse_argument.h>
#include <seqan/arg_parse/arg_parse_exceptions.h>

namespace seqan {

// ----------------------------------------------------------------------------
// Class ArgParseOption
// ----------------------------------------------------------------------------

/*!
 * @class ArgParseOption
 * @extends ArgParseArgument
 * @headerfile <seqan/arg_parse.h>
 * @brief Store information for a specific command line option.
 *
 * @signature class ArgParseOption : ArgParseArgument;
 *
 * You can use the shortcut @link ArgumentParser#addOption @endlink for adding options to @link ArgumentParser @endlink
 * object.
 */

/**
.Class.ArgParseOption
..base:Class.ArgParseArgument
..cat:Miscellaneous
..summary:Stores information for a specific command line option.
..signature:ArgParseOption
..remarks:A @Class.ArgParseOption@ object can be added to a @Class.ArgumentParser@ via @Function.ArgumentParser#addOption@.
..include:seqan/arg_parse.h
..see:Class.ArgParseArgument
..see:Class.ArgumentParser
*/

/*!
 * @fn ArgParseOption::ArgParseOption
 * @brief Constructor.
 *
 * @signature ArgParseOption::ArgParseOption(shortName, longName, helpText, argumentType[, argumentLabel[, isList[, numValues]]]);
 *
 * @param shortName     The short name of the argument.
 * @param longName      The long name of the argument (<tt>std::string</tt>).
 * @param helpText      The text to display as help (<tt>std::string</tt>).
 * @param argumentType  The type of the argument (@link ArgParseArgument::ArgumentType @endlink).
 * @param argumentLabel The label for the value to use in the help display, e.g. 'INT' in '--value INT'
 *                      (<tt>std::string</tt>).
 * @param isList        Flag for whether this option can be given multiple times (<tt>bool</tt>, <tt>true</tt> for
 *                      allowing multiple values).
 * @param numValues     Number of command line arguments that each option should bind (<tt>unsigned</tt>).
 *
 * @section Short and Long Option Names
 *
 * Short option names are passed to the program following one dash (e.g. '-i', '-e').  Long optiongs are passed
 * following two dashes (e.g. '--num-threads', '--error-rate').  The leading one or two dashes are not passed to the
 * constructor of ArgParseOption, however.
 *
 * @section Naming Conventions
 *
 * <ul>
 * <li>Keep the short name short, ideally to one character.</tt>
 * <li>Only use lower-case characters in the long name, use dashes to separate words, only use alphanumeric characters
 *     and dashes, no underline characters</li>.
 * </ul>
 */

/**
.Memfunc.ArgParseOption#ArgParseOption
..class:Class.ArgParseOption
..summary:Constructor
..signature:ArgParseOption(shortName, longName, helpText, argumentType[, argumentLabel[, isList]])
..param.shortName:A std::string containing the short-name option identifier (e.g. $"h"$ for the $-h/--help$ option).
Although not suggested the short-name can contain more than 1 character.
...remarks:Note that the leading "-" is not passed.
..param.longName:A std::string containing the long-name option identifier (e.g. $"help"$ for the $-h/--help$ option).
...remarks:Note that the leading "--" is not passed.
..param.helpText:A std::string containing the help text associated with this option.
..param.argument:A $ArgParseArgument::ArgumentType$ for the option (e.g., an integer argument).
...type:Class.ArgParseArgument
..param.argumentLabel:The label to use for the argument in the help text, e.g. $"NUMBER"$ for an integer. Optional.
...default:$""$
...type:nolink:$char const *$
..param.isList:Whether or not the argument allows multiple values.
...default:$false$
...type:nolink:$bool$
*/

class ArgParseOption :
    public ArgParseArgument
{
public:
    // ----------------------------------------------------------------------------
    // Members to specify the names of the ArgParseOption
    // ----------------------------------------------------------------------------
    std::string         shortName;     // short option name
    std::string         longName;      // long option name

    // ----------------------------------------------------------------------------
    // Members representing type, content and restrictions of the ArgParseOption
    // ----------------------------------------------------------------------------
    bool                _isFlag;      // true if this a bool option, that has no
                                      // argument we will internally represent it as a
                                      // string option set to either "true" or "false"
    bool                _isRequired;  // true if this ArgParseOption must be set
    bool                _isHidden;    // true if this ArgParseOption should not be
                                      // shown on the command line

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help,
                   ArgumentType argumentType,
                   std::string const & argumentLabel = "",
                   bool isListArgument = false,
                   unsigned numberOfValues = 1) :
        ArgParseArgument(argumentType, argumentLabel, isListArgument, numberOfValues),
        shortName(_shortName),
        longName(_longName),
        _isFlag(false),
        _isRequired(false),
        _isHidden(false)
    {
        _helpText = _help;
    }

    ArgParseOption(std::string const & _shortName,
                   std::string const & _longName,
                   std::string const & _help) :
        ArgParseArgument(ArgParseArgument::STRING),
        shortName(_shortName),
        longName(_longName),
        _isFlag(true),
        _isRequired(false),
        _isHidden(false)
    {
        defaultValue.push_back("false");
        setValidValues(*this, "true false");
        _helpText = _help;
    }

};

// ----------------------------------------------------------------------------
// Function isStringArgument()
// ----------------------------------------------------------------------------

inline bool isStringArgument(ArgParseOption const & me)
{
    return isStringArgument(static_cast<ArgParseArgument>(me)) && !me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isBooleanOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#isBooleanOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether an option is a switch.
 *
 * @signature bool isBooleanOption(option);
 *
 * @param option The ArgParseOption object to query.
 *
 * @return bool <tt>true</tt> if <tt>option</tt> is a switch and <tt>false</tt> otherwise.
 *
 * Flags are given without arguments, e.g. the <tt>-l</tt> flag in the Unix <tt>ls</tt> command.
 */

/**
.Function.ArgParseOption#isBooleanOption
..class:Class.ArgParseOption
..summary:Returns whether option is a switch.
..cat:Miscellaneous
..signature:isBooleanOption(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is a switch.
..include:seqan/arg_parse.h
*/

inline bool isBooleanOption(ArgParseOption const & me)
{
    return me._isFlag;
}

// ----------------------------------------------------------------------------
// Function isHidden()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#isHidden
 * @headerfile <seqan/arg_parse.h>
 * @brief Return whether an option is hidden on the help screen.
 *
 * @signature bool isHidden(option);
 *
 * @param option The ArgParseOption object to query.
 *
 * @return bool <tt>true</tt> if it is hidden, <tt>false</tt> otherwise.
 *
 * By default, options are not hidden.
 */

/**
.Function.ArgParseOption#isHidden
..class:Class.ArgParseOption
..summary:Returns whether option is hidden on the help screen. Default is false.
..cat:Miscellaneous
..signature:isHidden(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is hidden on the help screen.
..include:seqan/arg_parse.h
*/

inline bool isHidden(ArgParseOption const & me)
{
    return me._isHidden;
}

// ----------------------------------------------------------------------------
// Function hideOption()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#hideOption
 * @headerfile <seqan/arg_parse.h>
 * @brief Hides the ArgParseOption from the help screen (or shows it again).
 *
 * @signature void hideOption(option[, hide]);
 *
 * @param option The ArgParseOption object to set the hidden flag of.
 * @param hide   <tt>bool</tt> that indicates whether to hide the flag (default: <tt>true</tt>)
 */

/**
.Function.ArgParseOption#hideOption
..class:Class.ArgParseOption
..summary:Hides the ArgParseOption from the help screen.
..cat:Miscellaneous
..signature:hideOption(option [, hide])
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.hide:The new visibility of the option. Default is false.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void hideOption(ArgParseOption & me, bool hide = true)
{
    me._isHidden = hide;
}

// ----------------------------------------------------------------------------
// Function isRequired()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#isRequired
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns whether the option is mandatory.
 *
 * @signature bool isRequired(option);
 *
 * @param option The ArgParseOption object to query.
 * 
 * @return bool <tt>true</tt> if the option is mandatory and <tt>false</tt> if it not.
 *
 * By default, options are not mandatory.
 */

/**
.Function.ArgParseOption#isRequired
..class:Class.ArgParseOption
..summary:Returns whether the option is mandatory.
..cat:Miscellaneous
..signature:isRequired(option)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..returns:$true$ if the option is mandatory.
..include:seqan/arg_parse.h
*/

inline bool isRequired(ArgParseOption const & me)
{
    return me._isRequired;
}

// ----------------------------------------------------------------------------
// Function setDefaultValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#setDefaultValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Set the default value for the given option.
 *
 * @signature void setDefaultValue(option, v);
 *
 * @param option The ArgParseOption to set the default value for.
 * @param v      The value to set, (any type that can be streamed into an <tt>std::stringstream</tt>).
 */

// TODO(holtgrew): Deprecate in favour of string-only variant?

/**
.Function.ArgParseOption#setDefaultValue
..summary:Sets the default value for the given option.
..cat:Miscellaneous
..remarks:Note that this overwrites any previously given default values.
..signature:setDefaultValue(option, value)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.value:The new default value.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline void setDefaultValue(ArgParseOption & me, const TValue & value)
{
    try
    {
        std::stringstream strm;
        strm << value;

        // check if all constraints are satisfied
        _checkValue(me, strm.str());

        // clear old values
        me.defaultValue.clear();

        // add defaultValue
        me.defaultValue.push_back(strm.str());
    }
    catch (ParseException & ex)
    {
        SEQAN_FAIL("Default value does not satisfy the restrictions:\n %s", ex.what());
    }
}

// ----------------------------------------------------------------------------
// Function addDefaultValue()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#addDefaultValue
 * @headerfile <seqan/arg_parse.h>
 * @brief Adds/appends a new value to the list of default values.
 *
 * @signature void setDefaultValue(option, v);
 *
 * @param option The ArgParseOption to appen the default value for.
 * @param v      The value to append, (any type that can be streamed into an <tt>std::stringstream</tt>).
 * 
 * @section Remarks
 *
 * This function does not check any length restrictions for this value.
 */


/**
.Function.ArgParseOption#addDefaultValue
..summary:Adds/appends a new value to the list of default values.
..cat:Miscellaneous
..remarks:Note that this method does not check any length restrictions for this value.
..signature:addDefaultValue(option, value)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.value:The new default value.
..include:seqan/arg_parse.h
*/

template <typename TValue>
inline void addDefaultValue(ArgParseOption & me, const TValue & value)
{
    try
    {
        std::stringstream strm;
        strm << value;

        // check if all constraints are satisfied
        _checkValue(me, strm.str());

        // add defaultValue
        me.defaultValue.push_back(strm.str());
    }
    catch (ParseException & ex)
    {
        SEQAN_FAIL("Default value does not satisfy the restrictions:\n %s", ex.what());
    }
}

// ----------------------------------------------------------------------------
// Function setRequired()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#setRequired
 * @headerfile <seqan/arg_parse.h>
 * @brief Sets whether or not the option is mandatory.
 *
 * @signature void setRequired(option, required);
 *
 * @param option   The ArgParseOption to modify.
 * @param required Flag whether the option is mandatory (<tt>bool</tt>).
 *
 * By default, options are not mandatory.
 */

/**
.Function.ArgParseOption#setRequired
..class:Class.ArgParseOption
..summary:Sets whether or not the option is mandatory.
..cat:Miscellaneous
..signature:setRequired(option, required)
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..param.required:The new required value of the option.
...type:nolink:bool
..include:seqan/arg_parse.h
*/

inline void setRequired(ArgParseOption & me, bool required)
{
    me._isRequired = required;
}

// ----------------------------------------------------------------------------
// Function getArgumentLabel()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#getArgumentLabel
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns the argument label.
 *
 * @signature std::string getArgumentLabel(option);
 *
 * @param option The ArgParseOption object to query.
 *
 * @return std::string The argument label string.
 */

inline std::string const getArgumentLabel(ArgParseOption const & me)
{
    if (isBooleanOption(me))
        return "";
    else
        return getArgumentLabel(static_cast<ArgParseArgument>(me));
}

// ----------------------------------------------------------------------------
// Function getOptionName()                                    [ArgParseOption]
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#getOptionName
 * @headerfile <seqan/arg_parse.h>
 * @brief Returns the name of the ArgParseOption in a well-formatted way.
 *
 * @signature std::string getOptionName(option);
 *
 * @param option The ArgParseOption object to query.
 *
 * @return std::string The option name string.
 */

/**
 .Function.ArgParseOption#getOptionName
 ..class:Class.ArgParseOption
 ..summary:Returns the name of the @Class.ArgParseOption@ in a well formated way.
 ..cat:Miscellaneous
 ..signature:getOptionName(option)
 ..param.option:The @Class.ArgParseOption@ object.
 ...type:Class.ArgParseOption
 ..include:seqan/arg_parse.h
 ..returns:The name of the option as well formated string (e.g., -h, --help).
 */
inline std::string getOptionName(ArgParseOption const & me)
{
    std::stringstream stream;

    stream << (empty(me.shortName) ? "" : "-");
    stream << me.shortName;
    stream << ((empty(me.shortName) || empty(me.longName) ? "" : ", "));
    if (!empty(me.longName))
    {
        stream << "--";
        stream << me.longName;
    }
    return stream.str();
}

// ----------------------------------------------------------------------------
// Function write()                                            [ArgParseOption]
// ----------------------------------------------------------------------------

/*!
 * @fn ArgParseOption#write
 * @headerfile <seqan/arg_parse.h>
 * @brief Writes the basic informationa bout the ArgParseOption to the provided stream.
 *
 * @signature void write(stream, option);
 *
 * @param stream The @link StreamConcept stream @endlink to write to.
 * @param option The ArgParseOption object to write to <tt>stream</tt>.
 */

/**
.Function.ArgParseOption#write
..class:Class.ArgParseOption
..summary:Writes the basic information about the @Class.ArgParseOption@ to the provided stream.
..cat:Miscellaneous
..signature:write(stream, option)
..param.stream:The target stream.
..param.option:The @Class.ArgParseOption@ object.
...type:Class.ArgParseOption
..include:seqan/arg_parse.h
*/

template <typename TStream>
inline void write(TStream & target, ArgParseOption const & me)
{
    streamPut(target, '\t');
    streamPut(target, getOptionName(me));
    streamPut(target, '\t');
    streamPut(target, '\t');
    streamPut(target, me._helpText);
}

// ----------------------------------------------------------------------------
// operator<<()                                                [ArgParseOption]
// ----------------------------------------------------------------------------

// TODO(holtgrew): We need to work out a consistent scheme with operator<<().

template <typename TStream>
inline TStream & operator<<(TStream & target, ArgParseOption const & source)
{

    write(target, source);
    return target;
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_OPTION_H_
