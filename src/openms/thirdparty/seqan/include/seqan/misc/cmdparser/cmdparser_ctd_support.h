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
// Author: Björn Kahlert <Bjoern.Kahlert@fu-berlin.de>
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_CTD_SUPPORT_H_
#define CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_CTD_SUPPORT_H_


namespace seqan {

// ----------------------------------------------------------------------------
// Function _join()
// ----------------------------------------------------------------------------

/*
 * joins all elements of the the passed StringSet into a single CharString
 * the provided delimiter is used to separate the single entries in the
 * resulting CharString
 */
template <typename TValue>
inline CharString
_join(StringSet<TValue> const & v, CharString const & delimiter)
{
    typedef typename Iterator<StringSet<TValue> const, Rooted>::Type TStringSetIterator;

    std::stringstream joined;
    for (TStringSetIterator it = begin(v); it != end(v); goNext(it))
    {
        if (it != begin(v))
            joined << delimiter;
        joined << *it;
    }
    return CharString(joined.str());
}

// ----------------------------------------------------------------------------
// Function _xmlEscape()
// ----------------------------------------------------------------------------

/*
 * make sure that the text we put into the XML does not break the XML
 * candidates are
 *  " -> &quot;
 *  ' -> &apos;
 *  & -> &amp;
 *  < -> &lt;
 *  > -> &gt;
 */
template <typename TSequence>
inline TSequence _xmlEscape(TSequence const & original)
{
    TSequence escaped;
    for (typename Iterator<TSequence const, Rooted>::Type ch  = begin(original); ch != end(original); goNext(ch))
    {
        if (value(ch) == '"')
            append(escaped, "&quot;");
        else if (value(ch) == '\'')
            append(escaped, "&apos;");
        else if (value(ch) == '&')
            append(escaped, "&amp;");
        else if (value(ch) == '<')
            append(escaped, "&lt;");
        else if (value(ch) == '>')
            append(escaped, "&gt;");
        else
            append(escaped, *ch);
    }
    return escaped;
}

// ----------------------------------------------------------------------------
// Function _addMinMaxRestrictions()
// ----------------------------------------------------------------------------

inline void
_addMinMaxRestrictions(StringSet<CharString> & restrictions, CommandLineOption const & opt)
{
    CharString minMaxRestriction = "";
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

// ----------------------------------------------------------------------------
// Function _addValidValuesRestrictions()
// ----------------------------------------------------------------------------

inline void
_addValidValuesRestrictions(StringSet<CharString> & restrictions, CommandLineOption const & opt)
{
    if (length(opt.validValues) != 0)
    {
        for (Iterator<StringSet<CharString> const, Rooted>::Type valid = begin(opt.validValues); valid != end(opt.validValues); goNext(valid))
        {
            // for files we set *.(Name of the format)
            if (isOutputFile(opt) || isInputFile(opt))
            {
                CharString filetype = "*.";
                append(filetype, *valid);
                appendValue(restrictions, filetype);
            }
            else
            {
                appendValue(restrictions, *valid);
            }
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
_includeInCTD(CommandLineOption const & opt)
{
    return !(opt.shortName == "h" || opt.shortName == "V" || opt.longName == "write-ctd" || (opt.shortName == "" && opt.longName == ""));
}

// ----------------------------------------------------------------------------
// Function writeCTD()
// ----------------------------------------------------------------------------

/*
.Function.writeCTD
..summary:Exports the app's interface description to a .ctd file.
..cat:Miscellaneous
..signature:writeCTD(parser)
..param.parser:The @Class.CommandLineParser@ object.
...type:Class.CommandLineParser
..include:seqan/misc/misc_cmdparser.h
*/

inline void
writeCTD(CommandLineParser const & me)
{
    typedef Iterator<CommandLineParser::TOptionMap const>::Type TOptionMapIterator;
    TOptionMapIterator optionMapIterator;

    // create file [appname].ctd in working directory
    CharString ctdfilename;
    getOptionValueLong(me, "write-ctd", ctdfilename);

    std::ofstream ctdfile;
    ctdfile.open(toCString(ctdfilename));
    ctdfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    ctdfile << "<tool status=\"external\">\n";
    ctdfile << "\t<name>" << _xmlEscape(me._appName) << "</name>\n";
    ctdfile << "\t<version>" << _xmlEscape(me._versionText) << "</version>\n";
    ctdfile << "\t<description><![CDATA[" << _xmlEscape(me._appName) << ".]]></description>\n";
    ctdfile << "\t<manual><![CDATA[" << _xmlEscape(me._appName) << ".]]></manual>\n"; // TODO: as soon as we have a more sophisticated documentation embedded into the CmdParser, we should at this here
    ctdfile << "\t<docurl>Direct links in docs</docurl>\n";
    ctdfile << "\t<category>SeqAn - Sequence Analaysis</category>\n";
    ctdfile << "\t<mapping><![CDATA[\n";

    for (optionMapIterator = begin(me.optionMap); optionMapIterator != end(me.optionMap); optionMapIterator++)
    {

        CommandLineOption const & opt = *optionMapIterator;
        // filter help, version and ctd_export
        if (!_includeInCTD(opt))
            continue;

        CharString optionName = (opt.shortName != "" ? opt.shortName : opt.longName);
        CharString flagName = (opt.shortName != "" ? "-" : "--");
        append(flagName, optionName);

        ctdfile << "<mapparam CLISwitch=\"" << flagName << "\" name=\"" << _xmlEscape(me._appName) << "." << optionName << "\"/>\n";
    }

    ctdfile << "]]></mapping>\n";
    ctdfile << "\t<PARAMETERS version=\"1.3\" xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/Param_1_3.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">" << std::endl;
    ctdfile << "\t\t<NODE name=\"" << _xmlEscape(me._appName) << "\" description=\"???\">" << std::endl;

    for (optionMapIterator = begin(me.optionMap); optionMapIterator != end(me.optionMap); optionMapIterator++)
    {
        CommandLineOption const & opt = *optionMapIterator;

        // filter help, version and ctd_export
        if (!_includeInCTD(opt))
            continue;

        // prefer short name for options
        CharString optionName = (opt.shortName != "" ? opt.shortName : opt.longName);

        CharString type;

        if (isStringOption(opt))
            type = "string";
        else if (isIntOption(opt))
            type = "int";
        else if (isDoubleOption(opt))
            type = "double";

        // set up tags
        StringSet<CharString> tags;
        if (isInputFile(opt))
        {
            appendValue(tags, "input file");
        }
        if (isOutputFile(opt))
        {
            appendValue(tags, "output file");
        }
        if (isOptionMandatory(*optionMapIterator))
        {
            appendValue(tags, "required");
        }

        // set up restrictions
        StringSet<CharString> restrictions;
        _addValidValuesRestrictions(restrictions, opt);
        _addMinMaxRestrictions(restrictions, opt);


        ctdfile << "\t\t\t<ITEM " <<
        "name=\"" << _xmlEscape(optionName) << "\" " <<
        "value=\"" << _xmlEscape(opt.defaultValue) << "\" " <<
        "type=\"" << type << "\" " <<
        "description=\"" << _xmlEscape(opt.helpText) << "\" " <<
        "tags=\"" << _xmlEscape(_join(tags, ",")) << "\" " <<
        "restrictions=\"" << _xmlEscape(_join(restrictions, ",")) << "\"" <<
        "/>" << std::endl;
    }

    ctdfile << "\t\t</NODE>" << std::endl;
    ctdfile << "\t</PARAMETERS>" << std::endl;
    ctdfile << "</tool>" << std::endl;

    ctdfile.close();
}

} // namespace seqan

#endif // CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_CTD_SUPPORT_H_
