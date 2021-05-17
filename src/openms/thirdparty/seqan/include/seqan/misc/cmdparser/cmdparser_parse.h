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

#ifndef CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_PARSE_H_
#define CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_PARSE_H_

#include <seqan/misc/cmdparser/cmdoption.h>
#include <seqan/misc/cmdparser/cmdparser.h>
#include <seqan/misc/cmdparser/cmdparser_ctd_support.h>

namespace seqan {

// ----------------------------------------------------------------------------
// Function parse()
// ----------------------------------------------------------------------------

/*
 .Function.parse:
 ..summary:Parses the command line.
 ..cat:Miscellaneous
 ..signature:parse(parser, argc, argv[, errorStream])
 ..param.parser:The @Class.CommandLineParser@ object.
 ...type:Class.CommandLineParser
 ..param.argc:Count of the objects on the command line.
 ..param.argv:Array of the different command line arguments ($const char *argv[]$).
 ..param.errorStream:A stream where error messages are sent to.
 ..remarks:Must be called before retrieving options or arguments.
 ..returns:$true$ if all required arguments are set and parseable and neither the help nor version argument is set.
 ..include:seqan/misc/misc_cmdparser.h
 */

template <typename TErrorStream>
bool parse(CommandLineParser & me, int argc, const char * argv[], TErrorStream & estream)
{
    // if the appName wasn't set .. parse from command line
    if (empty(me._appName))
        me._appName = _parseAppName(argv[0]);

    for (int argument_index = 1; argument_index < argc; ++argument_index)
    {
        if (argv[argument_index][0] == '-')  // this is possibly an option value
        {
            CharString inParam = argv[argument_index];
            unsigned len = length(inParam);

            if (len == 1)
            {
                streamPut(estream, me._appName);
                streamPut(estream, ": invalid option '-'\n");
                return false;
            }
            else if (inParam[1] != '-') // maybe a combination of multiple bool opts
            {
                for (unsigned s = 1; s < len; ++s)
                {
                    unsigned e = len;
                    for (; s < e; --e)
                    {
                        if (hasOption(me, infix(inParam, s, e)))
                        {
                            CommandLineOption & opt = getOption(me, infix(inParam, s, e));
                            s = --e;
                            if (isBooleanOption(opt))
                                _assignOptionValue(me, opt, "true", estream);
                            else
                            {
                                int firstArgIndex = 0;

                                if (e < len - 1)
                                {
                                    // Try getting the first option argument from the remaining characters
                                    // of this program argument. Use-case: immediately adjacent option
                                    // values without separating space, as in `-x1` instead of `-x 1`.
                                    if (!_assignOptionValue(me, opt, suffix(inParam, e + 1), 0, estream))
                                        return false;

                                    firstArgIndex = 1;
                                    s = len - 1;
                                }

                                if (argument_index + opt.argumentsPerOption - firstArgIndex < argc)
                                {
                                    for (int t = firstArgIndex; t < opt.argumentsPerOption; ++t)
                                        if (!_assignOptionValue(me, opt, argv[++argument_index], t, estream))
                                            return false;
                                }
                                else  // no value available
                                {
                                    _reportMissingArgument(me, opt, estream);
                                    return false;
                                }
                            }
                        }
                    }
                    if (s == e)
                    {
                        CharString invalidOpt("-");
                        append(invalidOpt, suffix(inParam, s));

                        _reportInvalidOption(me, invalidOpt, estream);
                        return false;
                    }
                }
            }
            else if (inParam[1] == '-')  // this is a long option
            {
                unsigned t = 2;
                CharString longOpt, val;
                for (; t < len && inParam[t] != '='; ++t)
                    appendValue(longOpt, inParam[t], Generous());
                if (t < len) // this one is a --name=value option
                    val = suffix(inParam, t + 1);

                // We might already have a value
                if (hasOption(me, longOpt))
                {
                    CommandLineOption & opt = getOption(me, longOpt);

                    if (!empty(val))
                    {
                        if (opt.argumentsPerOption == 1)
                        {
                            if (!_assignOptionValue(me, opt, val, estream))
                                return false;
                        }
                        else
                        {
                            _reportMissingArgument(me, opt, estream);
                            return false;
                        }
                    }
                    else if (isBooleanOption(opt))
                    {
                        _assignOptionValue(me, opt, "true", estream);
                    }
                    else if (argument_index + opt.argumentsPerOption < argc)
                    {
                        for (int t = 0; t < opt.argumentsPerOption; ++t)
                            if (!_assignOptionValue(me, opt, argv[++argument_index], t, estream))
                                return false;
                    }
                    else  // no value available
                    {
                        _reportMissingArgument(me, opt, estream);
                        return false;
                    }
                }
                else
                {
                    CharString invalidOpt("--");
                    append(invalidOpt, longOpt);
                    _reportInvalidOption(me, invalidOpt, estream);
                    return false;
                }
            }
        }
        else  // this seems to be a normal argument
        {
            appendValue(me._arguments, argv[argument_index]);
        }
    }
    if (hasOption(me, "version") && isSet(me, "version"))
    {
        printVersion(me, estream);
        return false;
    }
    if (hasOption(me, "write-ctd") && isSet(me, "write-ctd"))
    {
        writeCTD(me);
        return false;
    }
    if (isSet(me, "help"))
    {
        printHelp(me, estream);
        return false;
    }
    if (argc == 1 && me._requiredArguments > 0)
    {
        // print short help and exit
        printShortHelp(me, estream);
        return false;
    }

    return _allMandatorySet(me) && (length(me._arguments) >= me._requiredArguments);
}

inline bool
parse(CommandLineParser & me, int argc, const char * argv[])
{
    return parse(me, argc, argv, std::cerr);
}

} // namespace seqan

#endif // CORE_INCLUDE_SEQAN_MISC_CMDPARSER_CMDPARSER_PARSE_H_
