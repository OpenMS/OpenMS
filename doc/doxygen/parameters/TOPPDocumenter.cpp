// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Mathias Walzer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <QtCore/QProcess>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace OpenMS;
using namespace Internal;

void convertINI2HTML(const Param& p, ostream& os)
{
  // the .css file is included via the Header.html (see doc/doxygen/common/Header.html)
  // TODO add some general description on how to handle subsections, what each column means, what the tags mean, etc.
  os << "<div class=\"ini_global\">\n";
  os << "<div class=\"legend\">\n";
  os << "<b>Legend:</b><br>\n";
  os << " <div class=\"item item_required\">required parameter</div>\n";
  os << " <div class=\"item item_advanced\">advanced parameter</div>\n";
  os << "</div>\n";

  Param::ParamIterator it = p.begin();
  String indentation = "  ";

  while (it != p.end())
  {
    string key = it.getName();

    //write opened/closed nodes
    const std::vector<Param::ParamIterator::TraceInfo>& trace = it.getTrace();
    for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
    {
      if (it2->opened) //opened node
      {
        String d = it2->description;
        d.substitute("\n", "<br>");
        os << indentation
           << R"(<div class="node"><span class="node_name">)"
           // TODO replace/remove weird "(TOPPAS) instance 1" nodes that only confuse people.
           << (String().fillLeft('+', (UInt) indentation.size() / 2) + it2->name)
           << "</span><span class=\"node_description\">"
           << (d)
           << "</span></div>"
           << "\n";
        indentation += "  ";
      }
      else //closed node
      {
        indentation.resize(indentation.size() - 2);
        //os << indentation << "</div>" << "\n";
      }
    }

    //write item
    String s_attr;
    String s_req;
    if (it->tags.find("advanced") != it->tags.end())
      s_attr += " item_advanced"; // optionally add advanced class
    if (it->tags.find("required") != it->tags.end())
      s_req += " item_required"; // optionally add required class
    ParamValue::ValueType value_type = it->value.valueType();
    //write opening tag
    os << indentation
       << "<div class=\"item"
       << s_attr
       << "\"><span class=\"item_name"
       << s_req
       << "\" style=\"padding-left:"
       << indentation.size() * 4
       << "px;\">"
       << (it->name)
       << "</span><span class=\"item_value\">"
       << it->value.toString()
       << "</span>"
       << "\n";

    //replace all critical characters in description
    String d = it->description;
    d.substitute("\n", "<br>");
    os << "<span class=\"item_description\">" << (d) << "</span>";

    //tags
    String list;
    for (auto tag_it = it->tags.begin(); tag_it != it->tags.end(); ++tag_it)
    {
      if (*tag_it == "advanced")
        continue; // do not list "advanced" or "required" (this is done by color coding)
      if (*tag_it == "required")
        continue;
      if (!list.empty())
        list += ", ";
      list += *tag_it;
    }
    os << "<span class=\"item_tags\">" << (list) << "</span>";

    //restrictions
    String restrictions = "";
    switch (value_type)
    {
    case ParamValue::INT_VALUE:
    case ParamValue::INT_LIST:
    {
      // TODO think about doing the same infinity replacement
      // for default values. A single ":" looks weird.
      bool min_set = (it->min_int != -numeric_limits<Int>::max());
      bool max_set = (it->max_int != numeric_limits<Int>::max());
      if (max_set || min_set)
      {
        if (min_set)
          restrictions += String(it->min_int);
        else
          restrictions += "-&#8734;"; // infinity symbol
        restrictions += ':';
        if (max_set)
          restrictions += String(it->max_int);
        else
          restrictions += "&#8734;";
      }
    }
    break;

    case ParamValue::DOUBLE_VALUE:
    case ParamValue::DOUBLE_LIST:
    {
      bool min_set = (it->min_float != -numeric_limits<double>::max());
      bool max_set = (it->max_float != numeric_limits<double>::max());
      if (max_set || min_set)
      {
        if (min_set)
          restrictions += String(it->min_float);
        else
          restrictions += "-&#8734;"; // infinity symbol
        restrictions += ':';
        if (max_set)
          restrictions += String(it->max_float);
        else
          restrictions += "&#8734;";
      }
    }
    break;

    case ParamValue::STRING_VALUE:
    case ParamValue::STRING_LIST:
      if (!it->valid_strings.empty())
      {
        // make sure browsers can word wrap with additional whitespace
        // TODO: If param name is *modification* just add a link to 
        //  a page with all modifications otherwise you get a HUGE list.
        //  Also think about a different separator, in case the restrictions have commas.
        restrictions.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ", ");
      }
      break;

    default:
      break;
    }
    if (restrictions.empty())
      restrictions = " "; // create content, such that the cell gets an underline

    os << "<span class=\"item_restrictions\">" << restrictions << "</span>";

    os << "</div>"; // end div item

    ++it;
  }

  os << "</div>\n"; // end global div
}

bool generate(const ToolListType& tools, const String& prefix, const String& binary_directory)
{
  bool errors_occured = false;
  for (ToolListType::const_iterator it = tools.begin(); it != tools.end(); ++it)
  {
    //start process
    QProcess process;
    process.setProcessChannelMode(QProcess::MergedChannels);
    QStringList env = QProcess::systemEnvironment();
    env << String("COLUMNS=110").toQString(); // Add an environment variable (used by each TOPP tool to determine width of help text (see TOPPBase))
    process.setEnvironment(env);

    String command = binary_directory + it->first;
#if defined(__APPLE__)
    if (it->first == "TOPPView" || it->first == "TOPPAS")
    {
      command = binary_directory + it->first + ".app/Contents/MacOS/" + it->first;
    }
#endif
#ifdef OPENMS_WINDOWSPLATFORM
    command += ".exe"; // otherwise File::exists() will fail
#endif

    ofstream f((String("output/") + prefix + it->first + ".cli").c_str());
    if (!File::exists(command))
    {
      stringstream ss;
      ss << "Errors occurred while generating the command line documentation for " << it->first << "!" << endl;
      ss << "Tool could not be found at '" << command << "'\n " << command << endl;
      f << ss.str();
      cerr << ss.str();
      errors_occured = true;
      f.close();
      continue;
    }
    else
    {
      ExternalProcess ep([&](const String& s) { f << s; }, 
                         [&](const String& s) { f << s; });
      String error_msg;
      if (ep.run(command.toQString(), QStringList() << "--help", "", false, error_msg, ExternalProcess::IO_MODE::READ_WRITE)
            != ExternalProcess::RETURNSTATE::SUCCESS)
      { // error while generation cli docu
        stringstream ss;
        ss << "Errors occurred while generating the command line documentation for " << it->first << "!" << endl;
        ss << "Output was: \n";
        ep.setCallbacks([&](const String& s) { ss << s; }, [&](const String& s) { ss << s; });
        ep.run(command.toQString(), QStringList() << "--help", "", false, error_msg, ExternalProcess::IO_MODE::READ_WRITE);
        ss << "\nCommand line was: \n " << command << endl;
        f << ss.str();
        cerr << ss.str();
        errors_occured = true;
        f.close();
        continue;
      }
    }
    f.close();

    //////
    // get the INI file and convert it into HTML
    //////
    if (it->first != "GenericWrapper" && // does not support -write_ini without a type
        it->first != "TOPPView" && // do not support -write_ini
        it->first != "TOPPAS")
    {
      String tmp_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_" + it->first + ".ini";
      const auto ini_command_args = QStringList() << "-write_ini" << tmp_file.toQString();
      
      ExternalProcess ep([&](const String& s) { f << s; }, [&](const String& s) { f << s; });
      String error_msg;
      if (ep.run(command.toQString(), ini_command_args, "", false, error_msg,
                 ExternalProcess::IO_MODE::READ_WRITE)
            != ExternalProcess::RETURNSTATE::SUCCESS
          || ! File::exists(tmp_file))
      { // error while generation cli docu
        std::cerr << "Errors occurred while writing ini file for " << it->first << "!" << std::endl;
        std::cerr << "Command line was: \n " << command << ini_command_args.join(" ").toStdString() << std::endl;
        errors_occured = true;
        continue;
      }

      // load content of written ini file
      Param p;
      ParamXMLFile pf;
      pf.load(tmp_file, p);
      File::remove(tmp_file);
      ofstream f_html((String("output/") + prefix + it->first + ".html").c_str());
      convertINI2HTML(p, f_html);
      f_html.close();
    }
  }
  return errors_occured;
}

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    cerr << "Please specify the path where the TOPP binaries are located." << endl;
    return EXIT_FAILURE;
  }

  String binary_directory = String(argv[1]).ensureLastChar('/');

  if (!File::exists(binary_directory))
  {
    cerr << "The given binary directory does not exist. Aborting." << endl;
    return EXIT_FAILURE;
  }

  //TOPP tools
  ToolListType topp_tools = ToolHandler::getTOPPToolList(true); // include GenericWrapper (can be called with --help without error, even though it has a type)
  topp_tools["TOPPView"] = Internal::ToolDescription(); // these two need to be excluded from writing an INI file later!
  topp_tools["TOPPAS"] = Internal::ToolDescription();

  bool errors_occured = generate(topp_tools, "TOPP_", binary_directory);

  if (errors_occured)
  {
    // errors occurred while generating the TOPP CLI docu .. tell the user
    cerr << "Errors occurred while generating the command line documentation for some of the TOPP tools." << endl;
    return EXIT_FAILURE;
  }
  else
  {
    return EXIT_SUCCESS;
  }
}
