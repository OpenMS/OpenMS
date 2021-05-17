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
// Author: Bjoern Kahlert <Bjoern.Kahlert@fu-berlin.de>
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_

#include <seqan/sequence.h>

#include <seqan/arg_parse/xml_support.h>
#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_doc.h>

#include <fstream>

namespace seqan {

// ----------------------------------------------------------------------------
// Function _toText()
// ----------------------------------------------------------------------------
// Removes formatting (\fI, \fB, and \fP).
template <typename TSequence>
TSequence _toText(TSequence const & input)
{
    TSequence buffer = xmlEscape(input);
    TSequence result;
    String<TSequence> openTags;

    typedef typename Iterator<TSequence const, Standard>::Type TIterator;
    TIterator endIt = end(input, Standard());
    for (TIterator it = begin(input, Standard()); it != endIt; goNext(it))
    {
        if (*it == '\\')
        {
            // Handle escape sequence, we interpret only "\-", "\fI", and "\fB".
            goNext(it);
            SEQAN_ASSERT_NOT(it == endIt);
            if (*it == '-')
            {
                appendValue(result, *it);
            }
            else if (*it == 'f')
            {
                goNext(it);
                SEQAN_ASSERT_NOT(it == endIt);
                if (*it == 'I')
                {
                    appendValue(openTags, "i");
                }
                else if (*it == 'B')
                {
                    appendValue(openTags, "b");
                }
                else if (*it == 'P')
                {
                    SEQAN_ASSERT_NOT(empty(openTags));
                    eraseBack(openTags);
                }
                else
                {
                    append(result, "\\f");
                    appendValue(result, *it);
                }
            }
            else
            {
                appendValue(result, '\\');
                appendValue(result, *it);
            }
        }
        else
        {
            appendValue(result, *it);
        }
    }

    return result;
}

// ----------------------------------------------------------------------------
// Function _join()
// ----------------------------------------------------------------------------

/**
 * joins all elements of the the passed StringSet into a single CharString
 * the provided delimiter is used to separate the single entries in the
 * resulting CharString
 */
template <typename TValue>
inline std::string
_join(std::vector<TValue> const & v, std::string const & delimiter)
{
    typedef typename std::vector<TValue>::const_iterator TStringSetIterator;

    std::stringstream joined;
    for (TStringSetIterator it = v.begin(); it != v.end(); ++it)
    {
        if (it != v.begin())
            joined << delimiter;
        joined << *it;
    }
    return joined.str();
}

// ----------------------------------------------------------------------------
// Function _getPrefixedOptionName()
// ----------------------------------------------------------------------------

inline std::string
_getPrefixedOptionName(ArgParseOption const & opt)
{
    std::string optName = "";
    if (!empty(opt.longName))
        optName = "--" + opt.longName;
    else
        optName = "-" + opt.shortName;

    return optName;
}

// ----------------------------------------------------------------------------
// Function _getOptionName()
// ----------------------------------------------------------------------------

inline std::string
_getOptionName(ArgParseOption const & opt)
{
    if (!empty(opt.longName))
        return opt.longName;
    else
        return opt.shortName;
}

// ----------------------------------------------------------------------------
// Function _getRestrictions()
// ----------------------------------------------------------------------------
inline void
_getRestrictions(std::vector<std::string> & restrictions, ArgParseArgument const & opt)
{
    // we only extract non-file restrictions
    if (isOutputFileArgument(opt) || isInputFileArgument(opt))
        return;

    if (length(opt.validValues) != 0)
    {
        for (std::vector<std::string>::const_iterator valid = opt.validValues.begin();
             valid != opt.validValues.end();
             ++valid)
        {
            appendValue(restrictions, *valid);
        }
    }
    else
    {
        std::string minMaxRestriction = "";
        if (opt.minValue != "")
        {
            append(minMaxRestriction, opt.minValue);
            append(minMaxRestriction, ":");
        }
        if (opt.maxValue != "")
        {
            if (minMaxRestriction == "")
                append(minMaxRestriction, ":");
            append(minMaxRestriction, opt.maxValue);
        }

        if (minMaxRestriction != "")
            appendValue(restrictions, minMaxRestriction);
    }
}

// ----------------------------------------------------------------------------
// Function _addValidValuesRestrictions()
// ----------------------------------------------------------------------------

inline void
_getSupportedFormats(std::vector<std::string> & supported_formats, ArgParseArgument const & opt)
{
    // we check only file arguments
    if (!(isOutputFileArgument(opt) || isInputFileArgument(opt)))
        return;

    if (length(opt.validValues) != 0)
    {
        std::string filetype;
        for (std::vector<std::string>::const_iterator valid = opt.validValues.begin();
             valid != opt.validValues.end();
             ++valid)
        {
            SEQAN_ASSERT_NOT(empty(*valid));

            filetype = "*";

            // ensure . as separator between * and file-extension
            if (value(*valid, 0) != '.')
                appendValue(filetype, '.');

            append(filetype, *valid);
            appendValue(supported_formats, filetype);
        }
    }
}

// ----------------------------------------------------------------------------
// Function _includeInCTD()
// ----------------------------------------------------------------------------

/*
 * returns true if this option should be included in the ctd
 */
inline bool
_includeInCTD(ArgParseOption const & opt)
{
    return !(opt.longName == "help" || opt.longName == "version" || opt.longName == "write-ctd" || opt.longName == "export-help" || (opt.shortName == "" && opt.longName == ""));
}

// ----------------------------------------------------------------------------
// Function _indent()
// ----------------------------------------------------------------------------

inline std::string _indent(const int currentIndent)
{
    std::string indent = "";
    for (int i = 0; i < currentIndent; ++i)
        indent += "\t";
    return indent;
}

// ----------------------------------------------------------------------------
// Function _writeCLIElement()
// ----------------------------------------------------------------------------
inline void _writeCLIElement(std::ostream & ctdfile, int currentIndent, std::string const & optionIdentifier, std::string const & ref_name, bool isList)
{
    ctdfile << _indent(currentIndent)
            << "<clielement optionIdentifier=\"" << optionIdentifier
            << "\" isList=\"" << (isList ? "true" : "false") << "\">\n";

    ctdfile << _indent(currentIndent + 1) << "<mapping referenceName=\"" << ref_name << "\" />\n";

    ctdfile << _indent(currentIndent) << "</clielement>\n";
}

// ----------------------------------------------------------------------------
// Function _getManual()
// ----------------------------------------------------------------------------
inline std::string _getManual(ArgumentParser const & me)
{
    std::stringstream manual;
    for (unsigned i = 0; i < me._description.size(); ++i)
    {
        manual << _toText(me._description[i]) << std::endl;
    }
    return manual.str();
}

// ----------------------------------------------------------------------------
// Function writeCTD()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#writeCTD
 * @headerfile <seqan/arg_parse.h>\
 * @brief Export the app's interface description to a .ctd file.
 *
 * @signature bool writeCTD(parser[, stream]);
 *
 * @param parser The ArgumentParser to write the CTD file for.
 * @param stream A <tt>std::ostream</tt> to write to.  If omitted an output file with the name form the "write-ctd"
 *               parameter of the parser is used.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

/**
.Function.writeCTD
..summary:Exports the app's interface description to a .ctd file.
..cat:Miscellaneous
..signature:writeCTD(parser [, ctdfile])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.ctdfile:The stream where the ctd file will be written to. If non is given the function writes it to the file given in the write-ctd parameter.
..param.parser:The @Class.ArgumentParser@ object.
..returns:$true$ if the ctd file could be created correctly, $false$ otherwise.
..include:seqan/arg_parse.h
*/

inline bool
writeCTD(ArgumentParser const & me, std::ostream & ctdfile)
{
    typedef ArgumentParser::TOptionMap::const_iterator   TOptionMapIterator;
    typedef ArgumentParser::TArgumentMapSize TArgumentMapSize;

    ctdfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ctdfile << "<tool>\n";

    int currentIndent = 1;

    std::string toolname(toCString(xmlEscape(getAppName(me))));

    // remove "_" in the tool name and make the following letter uppercase
    std::string class_name;
    bool upcase = true;
    for (unsigned i = 0; i < toolname.size(); ++i)
    {
        if (toolname[i] == '_')
        {
            upcase = true;
            continue;
        }
        class_name.push_back(toolname[i]);
        if (upcase)
            class_name[class_name.size() - 1] = toupper(toolname[i]);
        upcase = false;
    }

    ctdfile << _indent(currentIndent) << "<name>" << class_name << "</name>\n";
    ctdfile << _indent(currentIndent) << "<executableName>" << toolname << "</executableName>\n";
    ctdfile << _indent(currentIndent) << "<version>" << xmlEscape(getVersion(me)) << "</version>\n";
    ctdfile << _indent(currentIndent) << "<description>" << xmlEscape(getShortDescription(me)) << "</description>\n";
    ctdfile << _indent(currentIndent) << "<manual>" << xmlEscape(_getManual(me)) << "</manual>\n"; // TODO(aiche): as soon as we have a more sophisticated documentation embedded into the CmdParser, we should at this here
    ctdfile << _indent(currentIndent) << "<docurl>http://www.seqan.de</docurl>\n";
    ctdfile << _indent(currentIndent) << "<category>" << xmlEscape(getCategory(me)) << "</category>\n";
    ctdfile << _indent(currentIndent++) << "<cli>\n";

    // the unix way 1st the options
    for (TOptionMapIterator optionMapIterator = me.optionMap.begin();
         optionMapIterator != me.optionMap.end();
         ++optionMapIterator)
    {
        ArgParseOption const & opt = *optionMapIterator;
        std::string optionIdentifier = _getPrefixedOptionName(opt);
        std::string refName = toolname + "." + _getOptionName(opt);

        if (_includeInCTD(opt))
        {
            _writeCLIElement(ctdfile, currentIndent, optionIdentifier, refName, isListArgument(opt));
        }
    }

    // add a warning to the CTD that arguments are hard to interpret by the users
    if (me.argumentList.size() > 0)
    {
        ctdfile << _indent(currentIndent)
                << "<!-- Following clielements are arguments."
                << " You should consider providing a help text to ease understanding. -->\n";
    }
    // then the arguments
    for (TArgumentMapSize argIdx = 0; argIdx != me.argumentList.size(); ++argIdx)
    {
        // arguments do not have an option identifier
        std::string optionIdentifier = "";
        std::stringstream refName;
        refName << toolname << "." << "argument-" << argIdx;
        _writeCLIElement(ctdfile, currentIndent, optionIdentifier, refName.str(), isListArgument(me.argumentList[argIdx]));
    }

    ctdfile << _indent(--currentIndent) << "</cli>\n";
    ctdfile << _indent(currentIndent++) << "<PARAMETERS  version=\"1.4\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_4.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";
    ctdfile << _indent(currentIndent++) << "<NODE name=\"" << toolname << "\" description=\"" << xmlEscape(getShortDescription(me)) << "\">\n";

    for (TOptionMapIterator optionMapIterator = me.optionMap.begin();
         optionMapIterator != me.optionMap.end();
         ++optionMapIterator)
    {
        ArgParseOption const & opt = *optionMapIterator;

        // exclude help, version, etc.
        if (!_includeInCTD(opt))
            continue;

        // prefer short name for options
        std::string optionName = _getOptionName(opt);

        std::string type;

        if (isStringArgument(opt) || isBooleanOption(opt))
            type = "string";
        else if (isIntegerArgument(opt) || isInt64Argument(opt))
            type = "int";
        else if (isDoubleArgument(opt))
            type = "double";

        // set up tags
        std::vector<std::string> tags;
        if (isInputFileArgument(opt))
        {
            appendValue(tags, "input file");
        }
        if (isOutputFileArgument(opt))
        {
            appendValue(tags, "output file");
        }
        if (isRequired(opt))
        {
            appendValue(tags, "required");
        }
        if (isHidden(opt))
        {
            appendValue(tags, "advanced");
        }

        // set up restrictions
        std::vector<std::string> restrictions;
        _getRestrictions(restrictions, opt);

        // set up supported formats
        std::vector<std::string> supported_formats;
        _getSupportedFormats(supported_formats, opt);

        if (isListArgument(opt))
        {
            ctdfile << _indent(currentIndent)
                    << "<ITEMLIST " << "name=\"" << xmlEscape(optionName) << "\" "
                    << "type=\"" << type << "\" "
                    << "description=\"" << xmlEscape(_toText(opt._helpText)) << "\" ";

            if (!empty(tags))
                ctdfile << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" ";
            if (!empty(restrictions))
                ctdfile << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\" ";
            if (!empty(supported_formats))
                ctdfile << "supported_formats=\"" << xmlEscape(_join(supported_formats, ",")) << "\" ";

            ctdfile << ">\n";

            for (size_t i = 0; i < opt.defaultValue.size(); ++i)
            {
                ctdfile << _indent(currentIndent + 1) << "<LISTITEM value=\"" << xmlEscape(opt.defaultValue[i]) << "\"/>\n";
            }
            ctdfile << _indent(currentIndent) << "</ITEMLIST>\n";
        }
        else
        {
            ctdfile << _indent(currentIndent)
                    << "<ITEM " << "name=\"" << xmlEscape(optionName) << "\" "
                    << "value=\"" << xmlEscape(_join(opt.defaultValue, ",")) << "\" "
                    << "type=\"" << type << "\" "
                    << "description=\"" << xmlEscape(_toText(opt._helpText)) << "\" ";

            if (!empty(tags))
                ctdfile << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" ";
            if (!empty(restrictions))
                ctdfile << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\" ";
            if (!empty(supported_formats))
                ctdfile << "supported_formats=\"" << xmlEscape(_join(supported_formats, ",")) << "\" ";

            ctdfile << " />\n";
        }
    }

    for (TArgumentMapSize argIdx = 0; argIdx != me.argumentList.size(); ++argIdx)
    {
        ArgParseArgument arg = me.argumentList[argIdx];

        // prefer short name for options
        std::stringstream argumentNameStream;
        argumentNameStream << "argument-" << argIdx;
        std::string optionName = argumentNameStream.str();

        std::string type;

        if (isStringArgument(arg))
            type = "string";
        else if (isIntegerArgument(arg) || isInt64Argument(arg))
            type = "int";
        else if (isDoubleArgument(arg))
            type = "double";

        // set up tags
        std::vector<std::string> tags;
        appendValue(tags, "required");
        if (isInputFileArgument(arg))
        {
            appendValue(tags, "input file");
        }
        if (isOutputFileArgument(arg))
        {
            appendValue(tags, "output file");
        }

        // set up restrictions
        std::vector<std::string> restrictions;
        _getRestrictions(restrictions, arg);

        // set up supported formats
        std::vector<std::string> supported_formats;
        _getSupportedFormats(supported_formats, arg);

        ctdfile << _indent(currentIndent)
                << "<ITEM" << (isListArgument(arg) ? "LIST" : "") << " name=\"" << xmlEscape(optionName) << "\" "
                << (isListArgument(arg) ? " " : "value=\"\" ")
                << "type=\"" << type << "\" "
                << "description=\"" << xmlEscape(_toText(arg._helpText)) << "\" "; // it will be "" in most cases but we try
        if (!empty(tags))
            ctdfile << "tags=\"" << xmlEscape(_join(tags, ",")) << "\" ";
        if (!empty(restrictions))
            ctdfile << "restrictions=\"" << xmlEscape(_join(restrictions, ",")) << "\" ";
        if (!empty(supported_formats))
            ctdfile << "supported_formats=\"" << xmlEscape(_join(supported_formats, ",")) << "\" ";

        ctdfile << " />\n";
    }

    ctdfile << _indent(--currentIndent) << "</NODE>\n";
    ctdfile << _indent(--currentIndent) << "</PARAMETERS>\n";
    ctdfile << "</tool>" << std::endl;

    return true;
}

inline bool
writeCTD(ArgumentParser const & me)
{
    // create file [appname].ctd in working directory
    std::string ctdfilename;
    getOptionValue(ctdfilename, me, "write-ctd");

    std::ofstream ctdfile;
    ctdfile.open(toCString(ctdfilename));

    if (!ctdfile.is_open())
    {
        std::cerr << getAppName(me) << ": Unable to create ctd file: " << ctdfilename << std::endl;
        return false;
    }

    writeCTD(me, ctdfile);

    ctdfile.close();
    return true;
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_CTD_SUPPORT_H_
