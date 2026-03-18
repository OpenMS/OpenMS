// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>


#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QFileInfo>
#include <QDir>
#include <QRegularExpression>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_GenericWrapper GenericWrapper

@brief Allows generically the wrapping of external tools.
<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; GenericWrapper &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any file the external tool can read </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool reading the output format </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (to produce pepXML) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> &rarr; GenericWrapper (type 'ProteinProphet') &rarr;</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (protXML to idXML) </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> RAW file </td>
            <td VALIGN="middle" ROWSPAN=1> &rarr; GenericWrapper (type 'RAWFileConvert') &rarr;</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool accepting mzML </td>
        </tr>

    </table>
</CENTER>

This tool is a wrapper to call external (non-OpenMS) executables/scripts.
Each supported tool is represented by a certain <tt>type</tt>.
Each type exposes certain parameters which you can set (usually at least a <tt>in</tt> and <tt>out</tt>).

To obtain support for more external programs, visit the OpenMS website or (if you cannot find your tool there) ask on the OpenMS mailing list.

<b>The following section is for experts only, who want to add their own external tool:</b>

Each external tool is configured via a wrapper XML file in 'OpenMS/share/OpenMS/TOOLS/EXTERNAL'. All files have the ending .ttd (TOPP tool description).
You can add one or more wrappers (i.e. types) per file, but we recommend one. The filename does not really matter, but it should be descriptive.

The ttd file has the following structure:
<table>

<tr><th>type</th><td>
The name of the type which is added to list of valid GenericWrapper types. It should be unique, otherwise you get a fatal error.
</td></tr>

<tr><th>category</th><td>
Category for TOPPAS.
</td></tr>

<tr><th>cloptions</th><td>
Command line options (arguments) appended to the executable.
This string might contain placeholders of the form "%&lt;i&gt;"
where each placeholder will be substituted with a value that is determined in the
mappings section (see below).

Example:
@code
  <cloptions>-o "%1" --mzML "%2"</cloptions>
@endcode
</td></tr>

<tr><th>path</th><td>
Path (can be relative) to the executable that is executed.
</td></tr>

<tr><th>mappings</th><td>
Used to replace placeholders with input parameters.
The mapping id corresponds to the placeholder in <tt>cloptions</tt>.
The template used as starting string is given in <tt>cl</tt>.
All tokens therein will be replaced and the result will be patched into the <tt>cloptions</tt> string.
Allowed tokens are:
<ul>
<li>\%TMP  --> The current temp directory, fetched using File::getTempDirectory()
<li>\%DIR --> directory prefix, e.g.:, c:/tmp/mzfile.mzML gives 'c:/tmp'
<li>\%BASENAME[file] --> the basename of a file, e.g. c:/tmp/myfile.mzML gives 'myfile'
<li>\%RND --> generates a long random number, which can be used to generate unique directory or file names in a &lt;file_pre&gt; tag
<li>\%WORKINGDIR --> expands to the current working directory (default is '.'), settable by &lt;workingdirectory&gt; tag in the .ttd file.
<li>\%\%&lt;param&gt; --> any param registered in the ini_param section, e.g. '\%\%in'
</ul>

Example:
@code
  <mapping id="2" cl="-output_file %BASENAME[%%in].mgf -temp_dir %TMP -depth 3" />
@endcode
  </td></tr>

  <tr><th>ini_param</th><td>
  Contains part of a normal INI file with describes the parameters. Valid tags are those that are in the ParamXML scheme below 'NODE', e.g. 'ITEM'.
  Example:
@code
  <ITEM name="out" value="" type="string" description="output XML file containg regression line and confidence interval" tags="output file" />
  <ITEM name="mz_tolerance" value="1" type="float" description="Tolerance in m/z dimension" />
@endcode
</td></tr>

</table>

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_GenericWrapper.cli
<B>INI file documentation of this tool:</B>
*/
// no @htmlinclude TOPP_GenericWrapper.html since it needs a type to create an .INI (which would be only valid for this type)


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPGenericWrapper :
  public TOPPBase
{
public:
  TOPPGenericWrapper() :
    TOPPBase("GenericWrapper", "Allows the generic wrapping of external tools.")
  {
  }

protected:

  /**
    @brief format filenames and quote stringlists
  */
  String paramToString_(const Param::ParamEntry & p)
  {

    if (p.value.valueType() == ParamValue::STRING_LIST) // quote each element
    {
      StringList val = ListUtils::toStringList<std::string>(p.value);
      if (p.tags.count("input file") || p.tags.count("output file"))
      {
        for (Size i = 0; i < val.size(); ++i)
        {
          val[i] = QDir::toNativeSeparators(val[i].toQString());
        }
      }
      return "\"" + ListUtils::concatenate(val, "\" \"") + "\"";
    }
    if (p.tags.count("input file") || p.tags.count("output file"))
    {
      // ensure that file names are formated according to system spec
      return QDir::toNativeSeparators(String(p.value.toString()).toQString());
    }
    else
    {
      return p.value.toString();
    }
  }

  /**
    @brief Simple compare struct to sort a vector of String by the length of
    the contained strings

    */
  struct StringSizeLess
  {
    bool operator()(String const & left, String const & right) const
    {
      return left.size() < right.size();
    }

  };


  void createFragment_(String & fragment, const Param & param, const std::map<int, std::string>& optional_mappings = (std::map<int, std::string>()))
  {

    //std::cerr << "FRAGMENT: " << fragment << "\n\n";

    // e.g.:  -input %BASENAME[%%in].mzML

    // we have to make this little detour param -> vector<String>
    // to sort the param names by length, otherwise we have a
    // problem with parameter substitution
    // i.e., if A is a prefix of B and gets replaced first, the
    // suffix of B remains and will cause trouble, e.g.: "%%out" vs. "%%out_fm"
    vector<String> param_names;
    param_names.reserve(param.size());
    for (Param::ParamIterator it = param.begin(); it != param.end(); ++it)
    {
      param_names.push_back(it->name);
    }
    // sort by length
    std::sort(param_names.begin(), param_names.end(), [](auto &left, auto &right) {StringSizeLess cmp; return cmp(right, left);});

    // iterate through all input params and replace with values:
    SignedSize allowed_percent(0); // filenames might contain '%', which are allowed to remain there (and even must remain)
    for (vector<String>::iterator it = param_names.begin(); it != param_names.end(); ++it)
    {
      if (!fragment.hasSubstring("%%" + *it)) continue;

      String s_new = paramToString_(param.getEntry(*it));
      allowed_percent += s_new.length() - String(s_new).substitute("%", "").length();
      //std::cerr << "IN: " << s_new << "(" << allowed_percent << "\n";
      fragment.substitute("%%" + *it, s_new);
    }
    if (fragment.hasSubstring("%%"))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid '%%' found in '" + fragment + "' after replacing all parameters!", fragment);
    }
    // mapping replace> e.g.: %2
    // do it reverse, since %10 should precede %1
    for (std::map<int, std::string>::const_reverse_iterator it = optional_mappings.rbegin(); it != optional_mappings.rend(); ++it)
    {
      String m = String("%") + it->first;
      if (fragment.hasSubstring(m)) {
        writeDebug_(String("Replacing '") + m + "' in '" + fragment + "' by '" + it->second + "'\n", 10);
        fragment.substitute(m, it->second);
      }
    }

    // %TMP replace:
    fragment.substitute("%TMP", File::getTempDirectory());

    // %RND replace:
    fragment.substitute("%RND", String(UniqueIdGenerator::getUniqueId()));

    // %WORKINGDIR replace:
    fragment.substitute("%WORKINGDIR", tde_.working_directory);

    // %DIR% replace
    {
      QRegularExpression rx(R"(%DIR\[(.*)\])");
      rx.setPatternOptions(QRegularExpression::InvertedGreedinessOption);
      QString t_tmp = fragment.toQString();
      //std::cout << "fragment is:" << fragment << std::endl;
      for (const QRegularExpressionMatch& match : rx.globalMatch(fragment.toQString())) 
      {
        String value = match.captured(1);   // param name (hopefully)
        // replace in fragment:
        QFileInfo qfi(value.toQString());
        //std::cout << "match @ " << pos << " " << value << " --> " << qfi.canonicalPath() << "\n";
        t_tmp.replace(String("%DIR[" + value + "]").toQString(), qfi.canonicalPath());
      }
      fragment = t_tmp;
      //std::cout << "NEW fragment is:" << fragment << std::endl;
    }

    // %BASENAME% replace
    {
      QRegularExpression rx(R"(%BASENAME\[(.*)\])");
      rx.setPatternOptions(QRegularExpression::InvertedGreedinessOption);
      int count = 0;
      QString t_tmp = fragment.toQString();
      for (const QRegularExpressionMatch& match : rx.globalMatch(fragment.toQString())) 
      {
        //std::cout << "match @ " << pos << "\n";
        String value = match.captured(1); // param name (hopefully)
        // replace in fragment:
        QFileInfo qfi(value.toQString());
        //std::cout << "match @ " << pos << " " << value << " --> " << qfi.completeBaseName() << "\n";
        t_tmp.replace(String("%BASENAME[" + value + "]").toQString(), qfi.completeBaseName());
        ++count;
      }
      // update expected count of valid '%'
      allowed_percent -= (fragment.length() - String(fragment).substitute("%", "").length()) // original # of %
                         - (t_tmp.length() - String(t_tmp).substitute("%", "").length()) // new # of %
                         - count; // expected # of % due to %BASENAME
      fragment = String(t_tmp);
    }

    SignedSize diff = (fragment.length() - String(fragment).substitute("%", "").length()) - allowed_percent;
    //std::cerr << "allowed: " << allowed_percent << "\n" << "diff: " << diff << " in: " << fragment << "\n";
    if (diff > 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Mapping still contains '%' after substitution! Did you use % instead of %%?", fragment);
    }
    else if (diff < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Error: '%' from a filename where accidentally considered command tags! "
                                                                                              "This is a bug! Remove '%' from input filesnames to fix, but please report this as well!", fragment);
    }
    //std::cout << fragment << "'\n";
  }

  Internal::ToolExternalDetails tde_;

  ExitCodes wrapExit(const ExitCodes return_code) const
  {
    if (return_code != EXECUTION_OK)
    {
      OPENMS_LOG_ERROR << "\n" << tde_.text_fail << "\n";
    }
    return return_code;
  }

  void registerOptionsAndFlags_() override
  {
    registerSubsection_("ETool", "tool specific parameters");
    registerStringOption_("type", "", "", "Which external tool configuration to load?! See '" + ToolHandler::getExternalToolsPath() + "'.", true, false);
    setValidStrings_("type", ToolHandler::getTypes(toolName_()));
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    String type = getStringOption_("type"); // this will throw() if not set in param_
    // find params for 'type'
    Internal::ToolDescription gw = ToolHandler::getTOPPToolList(true)[toolName_()];
    for (Size i = 0; i < gw.types.size(); ++i)
    {
      if (type == gw.types[i])
      {
        return gw.external_details[i].param;
      }
    }
    // requested TDD is not found -- might be a custom TTD
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The value of 'Type' is invalid! Are you missing a TTD?", type);
  }

  ExitCodes main_(int, const char **) override
  {
    // find the config for the tool:
    String type = getStringOption_("type");


    Param tool_param = this->getParam_();

    // check required parameters (TOPPBase does not do this as we did not use registerInputFile_(...) etc)
    Param p = tool_param.copy("ETool:", true);
    for (Param::ParamIterator it = p.begin(); it != p.end(); ++it)
    {
      if ((it->tags).count("required") > 0)
      {
        String in = String(it->value.toString()).trim(); // will give '[]' for empty lists (hack, but DataValue class does not offer a convenient query)
        if (in.empty() || in == "[]") // any required parameter should have a value
        {
          OPENMS_LOG_ERROR << "The INI-parameter 'ETool:" << it->name << "' is required, but was not given! Aborting ..." << std::endl;
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
        else if ((it->tags).count("input file") > 0) // any required input file should exist
        {
          StringList ifs;
          switch (it->value.valueType())
          {
            case ParamValue::STRING_VALUE:
              ifs.push_back(it->value.toChar());
              break;
            case ParamValue::STRING_LIST:
              ifs = ListUtils::toStringList<std::string>(it->value);
              break;
            default:
              OPENMS_LOG_ERROR << "The INI-parameter 'ETool:" << it->name << "' is tagged as input file and thus must be a string! Aborting ...";
              return wrapExit(ILLEGAL_PARAMETERS);
          }
          for (StringList::const_iterator itf = ifs.begin(); itf != ifs.end(); ++itf)
          {
            if (!File::exists(*itf))
            {
              OPENMS_LOG_ERROR << "Input file '" << *itf << "' does not exist! Aborting ...";
              return wrapExit(INPUT_FILE_NOT_FOUND);
            }
          }
        }
      }
    }

    Internal::ToolDescription gw = ToolHandler::getTOPPToolList(true)[toolName_()];
    for (Size i = 0; i < gw.types.size(); ++i)
    {
      if (type == gw.types[i])
      {
        tde_ = gw.external_details[i];
        if (tde_.working_directory.trim().empty())
        {
          tde_.working_directory = ".";
        }
        break;
      }
    }

    OPENMS_LOG_INFO << tde_.text_startup << "\n";

    String command_args = tde_.commandline;
    // check for double spaces and warn
    if (command_args.hasSubstring("  "))
    {
      OPENMS_LOG_WARN << "Command line contains double spaces, which is not allowed. Condensing...\n";
      while (command_args.hasSubstring("  "))
      {
        command_args.substitute("  ", " ");
      }
      OPENMS_LOG_WARN << "result: " << command_args << std::endl;
    }

    writeDebug_("CommandLine from ttd (unprocessed): " + command_args, 1);

    // do "pre" moves (e.g. if the wrapped tool works on its data in-place (overwrites) it - we need to make a copy first
    // - we copy the file
    // - we set the value of the affected parameter to the copied tmp file, such that subsequent calls target the tmp file
    for (Size i = 0; i < tde_.tr_table.pre_moves.size(); ++i)
    {
      const Internal::FileMapping& fm = tde_.tr_table.pre_moves[i];
      // find target param:
      Param p = tool_param.copy("ETool:", true);
      String target = fm.target;
      if (!p.exists(target))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot find target parameter '" + target + "' being mapped from external tools output!", target);
      }
      String tmp_location = fm.location;
      // fragment's placeholder evaluation:

      createFragment_(tmp_location, p);

      // check if target already exists:
      String target_file = p.getValue(target).toString();
      if (File::exists(tmp_location))
      {
        if (!File::remove(tmp_location))
        {
          OPENMS_LOG_ERROR << "While writing a tmp file: Cannot remove conflicting file '" + tmp_location + "'. Check permissions! Aborting ...";
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
      }
      // create the temp file  tmp_location target_file
      writeDebug_(String("Copying '") + target_file + "' to '" + tmp_location + "'", 1);
      bool move_ok = QFile::copy(target_file.toQString(), tmp_location.toQString());
      if (!move_ok)
      {
        OPENMS_LOG_ERROR << "Copying the target file '" + tmp_location + "' from '" + target_file + "' failed! Aborting ...";
        return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
      }
      // set the input file's value to the temp file
      tool_param.setValue(String("ETool:") + target, tmp_location);
    }

    ///// construct the command line:
    std::map<int, std::string> mappings;  // remember the values for each mapping (for file_post substitution later on)
    // go through mappings (reverse because replacing %10 must come before %1):
    for (std::map<Int, String>::reverse_iterator it = tde_.tr_table.mapping.rbegin(); it != tde_.tr_table.mapping.rend(); ++it)
    {
      //std::cout << "mapping #" << it->first << "\n";
      String fragment = it->second;
      // fragment's placeholder evaluation:
      createFragment_(fragment, tool_param.copy("ETool:", true));

      // replace fragment in cl
      //std::cout << "replace : " << "%"+String(it->first) << " with '" << fragment << "\n";
      command_args.substitute("%" + String(it->first), fragment);

      // cache mapping
      mappings[it->first] = fragment;
    }

    QProcess builder;
    builder.setProcessChannelMode(QProcess::MergedChannels);
    String call = tde_.path + " " + command_args;

    writeDebug_("call command: " + call, 1);

    builder.setWorkingDirectory(tde_.working_directory.toQString());
    // TODO: start() with single argument is deprecated in Qt 5.15. Can probably be replaced with
    // QStringList commandArgs = QString::fromStdString(command_args).split(" ");
    // QString program = commandArgs.takeFirst();
    // builder.start(program, commandArgs);

    builder.start(call.toQString());

    if (!builder.waitForFinished(-1) || builder.exitStatus() != 0 || builder.exitCode() != 0)
    {
      OPENMS_LOG_ERROR << ("External tool returned with exit code (" + String(builder.exitCode()) + "), exit status (" + String(builder.exitStatus()) + ") or timed out. Aborting ...\n");
      OPENMS_LOG_ERROR << ("External tool output:\n" + String(QString(builder.readAll())));
      return wrapExit(EXTERNAL_PROGRAM_ERROR);
    }

    OPENMS_LOG_INFO << ("External tool output:\n" + String(QString(builder.readAll())));


    // post processing (file moving via 'file_post' command)
    for (Size i = 0; i < tde_.tr_table.post_moves.size(); ++i)
    {
      const Internal::FileMapping & fm = tde_.tr_table.post_moves[i];
      // find target param:
      Param p = tool_param.copy("ETool:", true);
      String source_file = fm.location;
      // fragment's placeholder evaluation:
      createFragment_(source_file, p, mappings);
      // check if target already exists:
      String target = fm.target;
      if (!p.exists(target))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot find target parameter '" + target + "' being mapped from external tools output!", target);
      }
      String target_file = p.getValue(target).toString();

      if (target_file.trim().empty())   // if target was not given, we skip the copying step (usually for optional parameters)
      {
        OPENMS_LOG_INFO << "Parameter '" + target + "' not given. Skipping forwarding of files.\n";
        continue;
      }
      // check if the target exists already (should not; if yes, delete it before overwriting it)
      if (File::exists(target_file))
      {
        if (!File::remove(target_file))
        {
          OPENMS_LOG_ERROR << "Cannot remove conflicting file '" + target_file + "'. Check permissions! Aborting ..." << std::endl;
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
      }
      // move to target
      writeDebug_(String("<file_post>: moving '") + source_file + "' to '" + target_file + "'", 1);
      if (!File::exists(source_file))
      {
        OPENMS_LOG_ERROR << "Moving the source file '" + source_file + "' during <file_post> failed, since it does not exist!\n"
                  << "Make sure the external program created the file and its filename is either\n"
                  << "unique or you only run one GenericWrapper at a time to avoid overwriting of files!\n"
                  << "Ideally, (if the external program allows to specify output filenames directly) avoid <file_post>\n"
                  << "in the TTD and request the output file directly. Aborting ..." << std::endl;
        return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
      }
      bool move_ok = QFile::rename(source_file.toQString(), target_file.toQString());
      if (!move_ok)
      {
        OPENMS_LOG_ERROR << "Moving the target file '" + target_file + "' from '" + source_file + "' failed!\n"
                  << "This file exists, but is either currently open for writing or otherwise blocked (concurrent process?). Aborting ..." << std::endl;
        return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
      }
    }

    OPENMS_LOG_INFO << tde_.text_finish << "\n";

    return wrapExit(EXECUTION_OK);
  }

};

int main(int argc, const char ** argv)
{
  TOPPGenericWrapper tool;
  return tool.main(argc, argv);
}

/// @endcond

