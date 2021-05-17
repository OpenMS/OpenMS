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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_PARSE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_PARSE_H_

#include <seqan/arg_parse/arg_parse_option.h>
#include <seqan/arg_parse/argument_parser.h>
#include <seqan/arg_parse/arg_parse_ctd_support.h>

namespace seqan {

// ----------------------------------------------------------------------------
// Function parse()
// ----------------------------------------------------------------------------

/*!
 * @fn ArgumentParser#parse
 * @headerfile <seqan/arg_parse.h>
 * @brief Parse command line parameters.
 *
 * @signature TResult parse(parser, argc, argv[, outStream, errStream]]);
 *
 * @param parser    The ArgumentParser to use for parsing and for storing parse results.
 * @param argc      The number of arguments (<tt>int</tt>).
 * @param argv      The arguments (<tt>const char * argv[]</tt>).
 * @param outStream The <tt>std::ostream</tt> to use for output.
 * @param errStream The <tt>std::ostream</tt> to use for error output.
 *
 * @return TResult The parse result, of type @link ArgumentParser::ParseResult @endlink.
 *
 * This function must be called before retrieving any options or arguments from the parser.
 */

/**
.Function.ArgumentParser#parse
..summary:Parses the command line.
..class:Class.ArgumentParser
..cat:Miscellaneous
..signature:parse(parser, argc, argv[, outputStream, errorStream])
..param.parser:The @Class.ArgumentParser@ object.
...type:Class.ArgumentParser
..param.argc:Count of the objects on the command line.
..param.argv:Array of the different command line arguments ($const char *argv[]$).
..param.errorStream:A stream where error messages are sent to.
..remarks:Must be called before retrieving options or arguments.
..returns:$true$ if all required arguments are set and parseable and neither the help nor version argument is set.
..include:seqan/arg_parse.h
*/

inline ArgumentParser::ParseResult parse(ArgumentParser & me,
                                         int argc,
                                         const char * argv[],
                                         std::ostream & outputStream,
                                         std::ostream & errorStream)
{
    typedef ArgumentParser::TArgumentMapSize TArgumentPosition;

    TArgumentPosition currentArgument = 0;

    // if the appName wasn't set .. parse from command line
    if (empty(getAppName(me)))
        _parseAppName(me, argv[0]);

    // we use exceptions here as an indicator for parse errors
    try
    {
        for (int arg = 1; arg < argc; ++arg)
        {
            if (argv[arg][0] == '-')  // this is possibly an option value
            {
                const std::string inParam = argv[arg];
                unsigned len = length(inParam);

                if (len == 1)
                {
                    throw InvalidOptionException("-");
                }
                else if (inParam[1] != '-') // maybe a combination of multiple bool opts
                {
                    for (unsigned s = 1; s < len; ++s)
                    {
                        unsigned e = len;
                        for (; s < e; --e)
                        {
                            if (hasOption(me, inParam.substr(s, e - s)))
                            {
                                ArgParseOption & opt = getOption(me, inParam.substr(s, e - s));
                                s = --e;
                                if (isBooleanOption(opt))
                                    _assignArgumentValue(opt, "true");
                                else
                                {
                                    if (e < len - 1)
                                    {
                                        std::stringstream what;
                                        what << "invalid combination of arguments -- " << inParam << std::endl;
                                        throw ParseException(what.str());
                                    }

                                    // assign the following values to this option
                                    if (arg + static_cast<int>(numberOfAllowedValues(opt)) < argc)
                                    {
                                        for (int t = 0; t < static_cast<int>(numberOfAllowedValues(opt)); ++t)
                                            _assignArgumentValue(opt, argv[++arg]);
                                    }
                                    else  // no value available
                                    {
                                        throw MissingArgumentException(opt.shortName);
                                    }
                                }
                            }
                        }
                        if (s == e)
                        {
                            throw InvalidOptionException(inParam.substr(s));
                        }
                    }
                }
                else if (inParam[1] == '-')  // this is a long option
                {
                    unsigned t = 2;
                    std::string longOpt, val;
                    for (; t < len && inParam[t] != '='; ++t)
                        longOpt += inParam[t];

                    if (t < len) // this one is a --name=value option
                        val = inParam.substr(t + 1);

                    // We might already have a value
                    if (hasOption(me, longOpt))
                    {
                        ArgParseOption & opt = getOption(me, longOpt);

                        if (!empty(val))
                        {
                            // we can only assign one value since it was set by --longOpt=val
                            if (numberOfAllowedValues(opt) == 1)
                                _assignArgumentValue(opt, val);
                            else
                                throw MissingArgumentException(longOpt);
                        }
                        else if (isBooleanOption(opt))
                            _assignArgumentValue(opt, "true");
                        else if (arg + static_cast<int>(numberOfAllowedValues(opt)) < argc)
                        {
                            for (int t = 0; t < static_cast<int>(numberOfAllowedValues(opt)); ++t)
                                _assignArgumentValue(opt, argv[++arg]);
                        }
                        else  // no value available
                        {
                            throw MissingArgumentException(longOpt);
                        }

                    }
                    else
                        throw InvalidOptionException(longOpt);
                }
            }
            else  // this seems to be a normal argument
            {
                // check if we have that much arguments
                if (me.argumentList.size() > currentArgument)
                {
                    ArgParseArgument & argument = getArgument(me, currentArgument);
                    _assignArgumentValue(argument, argv[arg]);

                    if (!isListArgument(argument))
                        ++currentArgument;
                }
                else
                {
                    throw ParseException("Too many arguments!");
                }
            }
        }
        if (hasOption(me, "version") && isSet(me, "version"))
        {
            printVersion(me, outputStream);
            return ArgumentParser::PARSE_VERSION;
        }
        if (hasOption(me, "write-ctd") && isSet(me, "write-ctd"))
        {
            if (writeCTD(me))
                return ArgumentParser::PARSE_WRITE_CTD;
            else
                return ArgumentParser::PARSE_ERROR;
        }
        if (isSet(me, "help"))
        {
            printHelp(me, outputStream);
            return ArgumentParser::PARSE_HELP;
        }
        if (isSet(me, "export-help"))
        {
            std::string format;
            getOptionValue(format, me, "export-help");
            printHelp(me, outputStream, format);
            return ArgumentParser::PARSE_EXPORT_HELP;
        }
        if (argc == 1 && (me.argumentList.size() > 0 || !_allRequiredSet(me)))
        {
            // print short help and exit
            printShortHelp(me, errorStream);
            return ArgumentParser::PARSE_HELP;
        }
    }
    catch (ParseException & ex)
    {
        errorStream << getAppName(me) << ": " << ex.what() << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    if (_allRequiredSet(me) && _allArgumentsSet(me))
        return ArgumentParser::PARSE_OK;
    else
    {
        // find missing options
        if (!_allRequiredSet(me))
        {
            for (unsigned o = 0; o < length(me.optionMap); ++o)
                if (!isSet(me.optionMap[o]) && isRequired(me.optionMap[o]))
                    errorStream << getAppName(me) << ": Missing value for option: " << getOptionName(me.optionMap[o]) << std::endl;
        }
        // and arguments
        if (!_allArgumentsSet(me))
        {
            errorStream << getAppName(me) << ": Not enough arguments were provided." << std::endl;
        }
        errorStream << "Try '" << getAppName(me) << " --help' for more information.\n";
        return ArgumentParser::PARSE_ERROR;
    }
}

inline ArgumentParser::ParseResult parse(ArgumentParser & me,
                                         int argc,
                                         const char * argv[])
{
    return parse(me, argc, argv, std::cout, std::cerr);
}

} // namespace seqan

#endif // SEQAN_CORE_INCLUDE_SEQAN_ARG_PARSE_ARG_PARSE_PARSE_H_
