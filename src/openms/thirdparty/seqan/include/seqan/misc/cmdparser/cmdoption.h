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

#ifndef CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDOPTION_H_
#define CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDOPTION_H_

#include <seqan/sequence.h>

namespace seqan {

/*
 * TODO: support some more formating options
 * TODO: store/return error code (invalid argument, invalid option, etc.)
 * TODO: support named arguments (e.g. <ARG1> -> <INPUT FILE>)
 * TODO: support user defined option argument names (e.g., -q QUALITY instead of -q NUM
 * TODO: correct parameter ordering of nearly all cmdparser function to be seqan conform f(out,in)
 */
struct OptionType
{
    // TODO(holtgrew): Should be all upper case!
    enum
    {
        Bool = 1,                   // option needs no argument, value is true iff given on command line
        Boolean = 1,                // option needs no argument, value is true iff given on command line
        String = 2,                 // argument is a string
        Int = 4,                    // ... an integer
        Integer = 4,                // ... an integer
        Double = 8,                 // ... a float
        Mandatory = 16,             // option must be set
        Label = 32,                 // automatically print a label for the argument(s) on the help screen
        List = 64,                  // option is a list of values
        Hidden = 128,               // hide this option from the help screen
        INPUTFILE = 256,            // this option is an input file .. is implicitly also a string, since paths/filenames are strings
        OUTPUTFILE = 512            // this option is an output file .. is implicitly also a string, since paths/filenames are strings
    };
};

// ----------------------------------------------------------------------------
// Class CommandLineOption
// ----------------------------------------------------------------------------

/*
.Class.CommandLineOption:
..cat:Miscellaneous
..summary:Stores information for a specific command line option.
..signature:CommandLineOption
..remarks:A @Class.CommandLineOption@ object can be added to a @Class.CommandLineParser@ via @Function.addOption@.
..include:seqan/misc/misc_cmdparser.h
*/

/*
.Memfunc.CommandLineOption#CommandLineOption:
..class:Class.CommandLineOption
..summary:Constructor
..signature:CommandLineOption ()
..signature:CommandLineOption (shortName, longName[, argumentsPerOption], helpText, optionType[, defaultValue])
..param.shortName:A @Shortcut.CharString@ containing the short-name option identifier (e.g. $"h"$ for the $-h/--help$ option).
Although not suggested the short-name can contain more than 1 character.
...remarks:Note that the leading "-" is not passed.
..param.longName:A @Shortcut.CharString@ containing the long-name option identifier (e.g. $"help"$ for the $-h/--help$ option).
...type:Shortcut.CharString
...remarks:Note that the leading "--" is not passed.
..param.argumentsPerOption:The number of required arguments per option (e.g. if set to 3 then 3 arguments must follow the option: "-foo x1 x2 x3").
...default:0 for boolean options and 1 for the rest.
..param.helpText:A @Shortcut.CharString@ containing the help text associated with this option.
...type:Shortcut.CharString
..param.optionType:Option type. This can be the sum of the some of the following values:
...tableheader:Flag|Value|Description
...table:$OptionType::Bool$ or $OptionType::Boolean$|1|Option needs no argument, value is true iff given on command line
...table:$OptionType::String$|2|Argument is a string
...table:$OptionType::Int$ or $OptionType::Integer$|4|An integer
...table:$OptionType::Double$|8|A float
...table:$OptionType::Mandatory$|16|Option must be set
...table:$OptionType::Label$|32|Automatically print a label for the argument(s) on the help screen
...table:$OptionType::List$|64|Option is a list of values
...table:$OptionType::Hidden$|128|Hide this option from the help screen
...table:$OptionType::INPUTFILE$|256|Argument is an input file
...table:$OptionType::OUTPUTFILE$|512|Argument is an output file
..param.defaultValue:The default value of this option.
...default:No default value.
*/

class CommandLineOption
{
public:
    CharString          longName;           // long option name
    CharString          shortName;          // short option name
    CharString          arguments;          // argument names seperated by spaces

    CharString          helpText;           // option description
    int                 optionType;         // option type
    int                 argumentsPerOption; // number of arguments per option

    // ----------------------------------------------------------------------------
    // Members to store the values
    // ----------------------------------------------------------------------------
    String<CharString>  defaultValue;
    String<CharString>  value;

    // ----------------------------------------------------------------------------
    // Members for restrictions
    // ----------------------------------------------------------------------------
    CharString            minValue;
    CharString            maxValue;
    StringSet<CharString> validValues;

    // ----------------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------------
    CommandLineOption() {}

    CommandLineOption(CharString const & _short,
                      CharString const & _long,
                      CharString const & _help,
                      int _type) :
        longName(_long),
        shortName(_short),
        helpText(_help),
        optionType(_type),
        argumentsPerOption(1),
        minValue(""),
        maxValue("")
    {}

    CommandLineOption(CharString const & _short,
                      CharString const & _long,
                      int _argumentsPerOption,
                      CharString const & _help,
                      int _type) :
        longName(_long),
        shortName(_short),
        helpText(_help),
        optionType(_type),
        argumentsPerOption(_argumentsPerOption),
        minValue(""),
        maxValue("")
    {}

    template <typename TValue>
    CommandLineOption(CharString const & _short,
                      CharString const & _long,
                      int _argumentsPerOption,
                      CharString const & _help,
                      int _type,
                      TValue const & _default) :
        longName(_long),
        shortName(_short),
        helpText(_help),
        optionType(_type),
        argumentsPerOption(_argumentsPerOption),
        minValue(""),
        maxValue("")
    {
        std::stringstream strm;
        strm << _default;
        appendValue(defaultValue, strm.str());
        append(helpText, " (default: ");
        append(helpText, strm.str());
        appendValue(helpText, ')');
    }

    template <typename TValue>
    CommandLineOption(CharString const & _short,
                      CharString const & _long,
                      CharString const & _help,
                      int _type,
                      TValue const & _default) :
        longName(_long),
        shortName(_short),
        helpText(_help),
        optionType(_type),
        argumentsPerOption(1),
        minValue(""),
        maxValue("")
    {
        std::stringstream strm;
        strm << _default;
        appendValue(defaultValue, strm.str());
        append(helpText, " (default: ");
        append(helpText, strm.str());
        appendValue(helpText, ')');
    }

};

// ----------------------------------------------------------------------------
// Function addArgumentText()
// ----------------------------------------------------------------------------

/*
.Function.addArgumentText:
..summary:Return a @Class.CommandLineOption@ object extended by an argument text.
..cat:Miscellaneous
..signature:addArgumentText(option, text)
..param.option:A @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.text:A @Shortcut.CharString@ containing the argument text.
...type:Shortcut.CharString
..returns:The option extended by the argument text.
Instead of using $option$, the return value can be used as argument for @Function.addOption@.
..remarks:The result type is a @Class.CommandLineOption@ object.
..include:seqan/misc/misc_cmdparser.h
*/

inline CommandLineOption
addArgumentText(CommandLineOption const & opt, CharString const & text)
{
    CommandLineOption temp = opt;
    temp.arguments = " ";
    append(temp.arguments, text);
    return temp;
}

// ----------------------------------------------------------------------------
// Function isStringOption()
// ----------------------------------------------------------------------------

/*
.Function.isStringOption
..summary:Returns whether option argument can be a string.
..cat:Miscellaneous
..signature:isStringOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option argument can be a string.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isStringOption(CommandLineOption const & me)
{
    return (me.optionType & (OptionType::String | OptionType::INPUTFILE | OptionType::OUTPUTFILE)) != 0;
}

// ----------------------------------------------------------------------------
// Function isBooleanOption()
// ----------------------------------------------------------------------------

/*
.Function.isBooleanOption
..summary:Returns whether option is a switch.
..cat:Miscellaneous
..signature:isBooleanOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option is a switch.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isBooleanOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Boolean) != 0;
}

// ----------------------------------------------------------------------------
// Function isDoubleOption()
// ----------------------------------------------------------------------------

/*
.Function.isDoubleOption
..summary:Returns whether option argument can be a double.
..cat:Miscellaneous
..signature:isDoubleOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option argument can be a double.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isDoubleOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Double) != 0;
}

// ----------------------------------------------------------------------------
// Function isIntOption()
// ----------------------------------------------------------------------------

/*
.Function.isIntOption
..summary:Returns whether option argument can be an integer.
..cat:Miscellaneous
..signature:isIntOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option argument can be an integer.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isIntOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Int) != 0;
}

// ----------------------------------------------------------------------------
// Function isHiddenOption()
// ----------------------------------------------------------------------------

/*
.Function.isHiddenOption
..summary:Returns whether option is hidden on the help screen.
..cat:Miscellaneous
..signature:isHiddenOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option is hidden on the help screen.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isHiddenOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Hidden) != 0;
}

// ----------------------------------------------------------------------------
// Function isOptionMandatory()
// ----------------------------------------------------------------------------

/*
.Function.isOptionMandatory
..summary:Returns whether option is mandatory.
..cat:Miscellaneous
..signature:isOptionMandatory(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option is mandatory.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isOptionMandatory(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Mandatory) != 0;
}

// ----------------------------------------------------------------------------
// Function isLabelOption()
// ----------------------------------------------------------------------------

/*
.Function.isLabelOption
..summary:Returns whether an option label should be printed on the help screen.
..cat:Miscellaneous
..signature:isLabelOption(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if an option label should be printed on the help screen.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isLabelOption(CommandLineOption const & me)
{
    return (me.optionType & OptionType::Label) != 0;
}

// ----------------------------------------------------------------------------
// Function isOptionList()
// ----------------------------------------------------------------------------

/*
.Function.isOptionList
..summary:Returns whether the option can be given multiple times.
..cat:Miscellaneous
..signature:isOptionList(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the option can be given multiple times on command line.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isOptionList(CommandLineOption const & me)
{
    return (me.optionType & OptionType::List) != 0;
}

// ----------------------------------------------------------------------------
// Function isInputFile()
// ----------------------------------------------------------------------------

/*
.Function.isInputFile
..summary:Returns whether the argument of the given option is an input file.
..cat:Miscellaneous
..signature:isInputFile(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the argument of the option is an input file.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isInputFile(CommandLineOption const & me)
{
    return (me.optionType & OptionType::INPUTFILE) != 0;
}

// ----------------------------------------------------------------------------
// Function isOutputFile()
// ----------------------------------------------------------------------------

/*
.Function.isOutputFile
..summary:Returns whether the argument of the given option is an output file.
..cat:Miscellaneous
..signature:isOutputFile(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:$true$ if the argument of the option is an output file.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline bool
isOutputFile(CommandLineOption const & me)
{
    return (me.optionType & OptionType::OUTPUTFILE) != 0;
}

// ----------------------------------------------------------------------------
// Function setOptionType()
// ----------------------------------------------------------------------------

/*
.Function.setOptionType:
..summary:Set the option type.
..cat:Miscellaneous
..signature:setOptionType(option, newOptionType)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.newOptionType:Option Type.
..see:Memfunc.CommandLineOption#CommandLineOption.param.optionType
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setOptionType(CommandLineOption & me, const int _newOptionType)
{
    me.optionType = _newOptionType;
}

// ----------------------------------------------------------------------------
// Function argumentText()
// ----------------------------------------------------------------------------

/*
.Function.argumentText
..summary:Returns the argument text of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:argumentText(option)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..returns:A text consisting of label and help text of the option.
...type:Shortcut.CharString
..include:seqan/misc/misc_cmdparser.h
*/

inline CharString
argumentText(CommandLineOption const & me)
{
    if (empty(me.arguments))
    {
        CharString label;
        if (isLabelOption(me))
        {
            if (isStringOption(me))
                label = " STR";
            else if (isIntOption(me) || isDoubleOption(me))
                label = " NUM";
            else if (isInputFile(me) || isOutputFile(me))
                label = " FILE";
            /*
            else if (isInputFileList(me) || isOutputFileList(me))
                label = " FILES";
            */
            if (me.argumentsPerOption >= 2)
            {
                std::stringstream strm;
                if (!empty(label))
                    for (int i = 0; i < me.argumentsPerOption; ++i)
                        strm << label << (i + 1);
                return strm.str();
            }
        }
        return label;
    }
    else
        return me.arguments;
}

// ----------------------------------------------------------------------------
// Helper Function _writeOptName()
// ----------------------------------------------------------------------------

template <typename TStream>
inline void
_writeOptName(TStream & target, CommandLineOption const & me)
{
    //IOREV _notio_ irrelevant for iorev
    streamPut(target, empty(me.shortName) ? "" : "-");
    streamPut(target, me.shortName);
    streamPut(target, (empty(me.shortName) || empty(me.longName)) ? "" : ", ");
    if (!empty(me.longName))
    {
        streamPut(target, "--");
        streamPut(target, me.longName);
    }
}

// ----------------------------------------------------------------------------
// Function write()                                           CommandLineOption
// ----------------------------------------------------------------------------

/*
.Function.write
..summary:Writes the basic information about the @Class.CommandLineOption@ to the provided stream.
..cat:Miscellaneous
..signature:write(stream,option)
..param.stream:The target stream.
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..include:seqan/misc/misc_cmdparser.h
*/

template <typename TStream>
inline void
write(TStream & target, CommandLineOption const & me)
{
    //IOREV _nodoc_ this specialization is not documented
    streamPut(target, '\t');
    _writeOptName(target, me);
    streamPut(target, '\t');
    streamPut(target, '\t');
    streamPut(target, me.helpText);
}

// ----------------------------------------------------------------------------
// operator<<()                                               CommandLineOption
// ----------------------------------------------------------------------------

template <typename TStream>
inline TStream &
operator<<(TStream & target, CommandLineOption const & source)
{
    //IOREV _nodoc_ this specialization is not documented
    write(target, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function setMinValue()
// ----------------------------------------------------------------------------

/*
.Function.setMinValue
..summary:Sets the minimum value of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:setMinValue(option,minValue)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.minValue:A @Shortcut.CharString@ containing a string representation of the minimum value of the @Class.CommandLineOption@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMinValue(CommandLineOption & me, const CharString _minValue)
{
    // TODO: check if value is applicable
    // TODO: if max exists, check if the interval is not empty
    me.minValue = _minValue;
}

// ----------------------------------------------------------------------------
// Function setMaxValue()
// ----------------------------------------------------------------------------

/*
.Function.setMaxValue
..summary:Sets the maximum value of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:setMaxValue(option,maxValue)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.maxValue:A @Shortcut.CharString@ containing a string representation of the maximum value of the @Class.CommandLineOption@.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setMaxValue(CommandLineOption & me, const CharString _maxValue)
{
    // TODO: check if value is applicable
    // TODO: if min exists, check if the interval is not empty
    me.maxValue = _maxValue;
}

// ----------------------------------------------------------------------------
// Function setValidValues()
// ----------------------------------------------------------------------------

/*
.Function.setValidValues
..summary:Sets the set of allowed values of a @Class.CommandLineOption@ object.
..cat:Miscellaneous
..signature:setValidValues(option,values)
..param.option:The @Class.CommandLineOption@ object.
...type:Class.CommandLineOption
..param.values:A $String<CharString>$ containing all valid entries for the option.
..include:seqan/misc/misc_cmdparser.h
*/

inline void
setValidValues(CommandLineOption & me, StringSet<CharString> const & _values)
{
    // TODO: check for collision with min/max
    // TODO: check if _values are applicable to option
    me.validValues = _values;
}

} // namespace seqan

#endif // CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDOPTION_H_
