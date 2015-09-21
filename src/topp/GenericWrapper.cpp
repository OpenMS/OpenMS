// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtCore/QProcess>
#include <QFileInfo>
#include <QDir>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
    @page TOPP_GenericWrapper GenericWrapper

    @brief Allows generically the wrapping of external tools.
<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ GenericWrapper \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any file the external tool can read </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool reading the output format </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (to produce pepXML) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> \f$ \longrightarrow \f$ GenericWrapper (type 'ProteinProphet') \f$ \longrightarrow \f$</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter (protXML to idXML) </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> RAW file </td>
            <td VALIGN="middle" ROWSPAN=1> \f$ \longrightarrow \f$ GenericWrapper (type 'RAWFileConvert') \f$ \longrightarrow \f$</td>
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
    @htmlinclude TOPP_GenericWrapper.html
*/


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

    if (p.value.valueType() == DataValue::STRING_LIST) // quote each element
    {
      StringList val = p.value;
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
      return QDir::toNativeSeparators(p.value.toQString());
    }
    else
    {
      return p.value;
    }
  }

  /**
    @brief Simple compare struct to sort a vector of String by the length of
    the contained strings

    */
  struct StringSizeLess :
    std::binary_function<String, String, bool>
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
    std::sort(param_names.begin(), param_names.end(), reverseComparator(StringSizeLess()));

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
    if (fragment.hasSubstring("%%")) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Invalid '%%' found in '" + fragment + "' after replacing all parameters!", fragment);

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
      QRegExp rx("%DIR\\[(.*)\\]");
      rx.setMinimal(true);
      int pos = 0;
      QString t_tmp = fragment.toQString();
      //std::cout << "fragment is:" << fragment << std::endl;
      while ((pos = rx.indexIn(t_tmp, pos)) != -1)
      {
        String value = rx.cap(1);   // param name (hopefully)
        // replace in fragment:
        QFileInfo qfi(value.toQString());
        //std::cout << "match @ " << pos << " " << value << " --> " << qfi.canonicalPath() << "\n";
        t_tmp = t_tmp.replace(String("%DIR[" + value + "]").toQString(), qfi.canonicalPath());
      }
      fragment = String(t_tmp);
      //std::cout << "NEW fragment is:" << fragment << std::endl;
    }

    // %BASENAME% replace
    {
      QRegExp rx("%BASENAME\\[(.*)\\]");
      rx.setMinimal(true);
      int pos = 0, count = 0;
      QString t_tmp = fragment.toQString();
      while ((pos = rx.indexIn(t_tmp, pos)) != -1)
      {
        //std::cout << "match @ " << pos << "\n";
        String value = rx.cap(1);   // param name (hopefully)
        // replace in fragment:
        QFileInfo qfi(value.toQString());
        //std::cout << "match @ " << pos << " " << value << " --> " << qfi.completeBaseName() << "\n";
        t_tmp = t_tmp.replace(String("%BASENAME[" + value + "]").toQString(), qfi.completeBaseName());
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
    if (diff > 0) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Mapping still contains '%' after substitution! Did you use % instead of %%?", fragment);
    else if (diff < 0) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Error: '%' from a filename where accidentally considered command tags! "
                                                                                              "This is a bug! Remove '%' from input filesnames to fix, but please report this as well!", fragment);

    //std::cout << fragment << "'\n";
  }

  Internal::ToolExternalDetails tde_;

  ExitCodes wrapExit(const ExitCodes return_code)
  {
    if (return_code != EXECUTION_OK)
    {
      LOG_ERROR << "\n" << tde_.text_fail << "\n";
    }
    return return_code;
  }

  void registerOptionsAndFlags_()
  {
    registerSubsection_("ETool", "tool specific parameters");
    registerStringOption_("type", "", "", "Which external tool configuration to load?! See '" + ToolHandler::getExternalToolsPath() + "'.", true, false);
    setValidStrings_("type", ToolHandler::getTypes(toolName_()));
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    String type = getStringOption_("type");
    // find params for 'type'
    Internal::ToolDescription gw = ToolHandler::getTOPPToolList(true)[toolName_()];
    for (Size i = 0; i < gw.types.size(); ++i)
    {
      if (type == gw.types[i])
      {
        return gw.external_details[i].param;
      }
    }
    return Param();
  }

  ExitCodes main_(int, const char **)
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
        if (it->value.toString().trim().empty())    // any required parameter should have a value
        {
          LOG_ERROR << "The INI-parameter '" + it->name + "' is required, but was not given! Aborting ...";
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
        else if ((it->tags).count("input file") > 0) // any required input file should exist
        {
          if (!File::exists(it->value))
          {
            LOG_ERROR << "Input file '" + String(it->value) + "' does not exist! Aborting ...";
            return wrapExit(INPUT_FILE_NOT_FOUND);
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
        if (tde_.working_directory.trim() == "") tde_.working_directory = ".";
        break;
      }
    }

    LOG_INFO << tde_.text_startup << "\n";

    String command_args = tde_.commandline;
    // check for double spaces and warn
    if (command_args.hasSubstring("  "))
    {
      LOG_WARN << "Commandline contains double spaces, which is not allowed. Condensing...\n";
      while (command_args.hasSubstring("  "))
      {
        command_args.substitute("  ", " ");
      }
      LOG_WARN << "result: " << command_args << std::endl;
    }

    writeDebug_("CommandLine from ttd (unprocessed): " + command_args, 1);

    // do "pre" moves (e.g. if the wrapped tool works on its data in-place (overwrites) it - we need to make a copy first
    // - we copy the file
    // - we set the value of the affected parameter to the copied tmp file, such that subsequent calls target the tmp file
    for (Size i = 0; i < tde_.tr_table.pre_moves.size(); ++i)
    {
      const Internal::FileMapping & fm = tde_.tr_table.pre_moves[i];
      // find target param:
      Param p = tool_param.copy("ETool:", true);
      String target = fm.target;
      if (!p.exists(target)) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot find target parameter '" + target + "' being mapped from external tools output!", target);
      String tmp_location = fm.location;
      // fragment's placeholder evaluation:

      createFragment_(tmp_location, p);

      // check if target already exists:
      String target_file = (String)p.getValue(target);
      if (File::exists(tmp_location))
      {
        if (!File::remove(tmp_location))
        {
          LOG_ERROR << "While writing a tmp file: Cannot remove conflicting file '" + tmp_location + "'. Check permissions! Aborting ...";
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
      }
      // create the temp file  tmp_location target_file
      writeDebug_(String("Copying '") + target_file + "' to '" + tmp_location + "'", 1);
      bool move_ok = QFile::copy(target_file.toQString(), tmp_location.toQString());
      if (!move_ok)
      {
        LOG_ERROR << "Copying the target file '" + tmp_location + "' from '" + target_file + "' failed! Aborting ...";
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
    builder.start(call.toQString());

    if (!builder.waitForFinished(-1) || builder.exitStatus() != 0 || builder.exitCode() != 0)
    {
      LOG_ERROR << ("External tool returned with non-zero exit code (" + String(builder.exitCode()) + "), exit status (" + String(builder.exitStatus()) + ") or timed out. Aborting ...\n");
      LOG_ERROR << ("External tool output:\n" + String(QString(builder.readAll())));
      return wrapExit(EXTERNAL_PROGRAM_ERROR);
    }

    LOG_INFO << ("External tool output:\n" + String(QString(builder.readAll())));


    // post processing (file moving via 'file_post' command)
    for (Size i = 0; i < tde_.tr_table.post_moves.size(); ++i)
    {
      const Internal::FileMapping & fm = tde_.tr_table.post_moves[i];
      // find target param:
      Param p = tool_param.copy("ETool:", true);
      String source = fm.location;
      // fragment's placeholder evaluation:
      createFragment_(source, p, mappings);
      // check if target already exists:
      String target = fm.target;
      if (!p.exists(target)) throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Cannot find target parameter '" + target + "' being mapped from external tools output!", target);
      String target_file = (String)p.getValue(target);

      if (target_file.trim().empty())   // if target was not given, we skip the copying step (usually for optional parameters)
      {
        LOG_INFO << "Parameter '" + target + "' not given. Skipping forwarding of files.\n";
        continue;
      }
      // check if the target exists already (should not; if yes, delete it before overwriting it)
      if (File::exists(target_file))
      {
        if (!File::remove(target_file))
        {
          LOG_ERROR << "Cannot remove conflicting file '" + target_file + "'. Check permissions! Aborting ..." << std::endl;
          return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
        }
      }
      // move to target
      writeDebug_(String("moving '") + source + "' to '" + target_file + "'", 1);
      bool move_ok = QFile::rename(source.toQString(), target_file.toQString());
      if (!move_ok)
      {
        LOG_ERROR << "Moving the target file '" + target_file + "' from '" + source + "' failed! Aborting ..." << std::endl;
        return wrapExit(CANNOT_WRITE_OUTPUT_FILE);
      }
    }

    LOG_INFO << tde_.text_finish << "\n";

    return wrapExit(EXECUTION_OK);
  }

};

int main(int argc, const char ** argv)
{
  TOPPGenericWrapper tool;
  return tool.main(argc, argv);
}

/// @endcond

