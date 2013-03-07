// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm, Mathias Walzer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/ToolHandler.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QtCore/QProcess>
#include <QDir>

#include <fstream>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace OpenMS;
using namespace Internal;

void convertINI2HTML(const Param & p, ostream & os)
{
  // the .css file is included via the Header.html (see doc/doxygen/common/Header.html)
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
    const std::vector<Param::ParamIterator::TraceInfo> & trace = it.getTrace();
    for (std::vector<Param::ParamIterator::TraceInfo>::const_iterator it2 = trace.begin(); it2 != trace.end(); ++it2)
    {
      if (it2->opened) //opened node
      {
        String d = it2->description;
        d.substitute("\n", "<br>");
        os << indentation  << "<div class=\"node\"><span class=\"node_name\">" << (String().fillLeft('+', (UInt) indentation.size() / 2) + it2->name) << "</span><span class=\"node_description\">" << (d) << "</span></div>" << "\n";
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
      s_attr += " item_advanced";                                                // optionally add advanced class
    if (it->tags.find("required") != it->tags.end())
      s_req += " item_required";                                                // optionally add required class
    DataValue::DataType value_type = it->value.valueType();
    //write opening tag
    os << indentation << "<div class=\"item" + s_attr + "\"><span class=\"item_name" + s_req + "\" style=\"padding-left:" << indentation.size() * 4 << "px;\">" << (it->name) << "</span><span class=\"item_value\">" << it->value.toString() << "</span>" << "\n";

    //replace all critical characters in description
    String d = it->description;
    d.substitute("\n", "<br>");
    os << "<span class=\"item_description\">" << (d) << "</span>";

    //tags
    String list;
    for (set<String>::const_iterator tag_it = it->tags.begin(); tag_it != it->tags.end(); ++tag_it)
    {
      if (*tag_it == "advanced")
        continue;                        // do not list "advanced" or "required" (this is done by color coding)
      if (*tag_it == "required")
        continue;
      if (!list.empty())
        list += ",";
      list += *tag_it;
    }
    os << "<span class=\"item_tags\">" << (list) << "</span>";

    //restrictions
    String restrictions = "";
    bool escape_restrictions(true);
    switch (value_type)
    {
    case DataValue::INT_VALUE:
    case DataValue::INT_LIST:
    {
      bool min_set = (it->min_int != -numeric_limits<Int>::max());
      bool max_set = (it->max_int != numeric_limits<Int>::max());
      if (max_set || min_set)
      {
        if (min_set)
          restrictions += String(it->min_int);
        else
          restrictions += "-&#8734;";        // infinity symbol
        restrictions += ':';
        if (max_set)
          restrictions += String(it->max_int);
        else
          restrictions += "&#8734;";
      }
      escape_restrictions = false;     // prevent html escape of infinity symbol
    }
    break;

    case DataValue::DOUBLE_VALUE:
    case DataValue::DOUBLE_LIST:
    {
      bool min_set = (it->min_float != -numeric_limits<DoubleReal>::max());
      bool max_set = (it->max_float != numeric_limits<DoubleReal>::max());
      if (max_set || min_set)
      {
        if (min_set)
          restrictions += String(it->min_float);
        else
          restrictions += "-&#8734;";        // infinity symbol
        restrictions += ':';
        if (max_set)
          restrictions += String(it->max_float);
        else
          restrictions += "&#8734;";
      }
      escape_restrictions = false;     // prevent html escape of infinity symbol
    }
    break;

    case DataValue::STRING_VALUE:
    case DataValue::STRING_LIST:
      if (it->valid_strings.size() != 0)
      {
        restrictions.concatenate(it->valid_strings.begin(), it->valid_strings.end(), ",");
      }
      break;

    default:
      break;
    }
    if (restrictions.empty())
      restrictions = " ";                       // create content, such that the cell gets an underline

    os << "<span class=\"item_restrictions\">" << restrictions << "</span>";

    os << "</div>"; // end div item

    ++it;
  }

  os << "</div>\n";  // end global div
}

bool generate(const ToolListType & tools, const String & prefix, const String & binary_directory)
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
    process.start(String(command + " --help").toQString());
    process.waitForFinished();

    ofstream f((String("output/") + prefix + it->first + ".cli").c_str());
    std::string lines = QString(process.readAll()).toStdString();
    if (process.error() != QProcess::UnknownError)
    {
      // error while generation cli docu
      f << "Errors occurred while generating the command line documentation for " << it->first << "!" << endl;
      f << "Output was: \n" << lines << endl;
      f << "Command line was: \n " << command << endl;
      errors_occured = true;
    }
    else
    {
      // write output
      f << lines;
    }
    f.close();

    //////
    // get the INI file and convert it into HTML
    //////
    if (it->first == "GenericWrapper")
      continue;                                  // does not support -write_ini without a type
    if (it->first == "TOPPView")
      continue;                            // does not support -write_ini
    if (it->first == "TOPPAS")
      continue;                          // does not support -write_ini
    String tmp_file = File::getTempDirectory() + "/" + File::getUniqueName() + "_" + it->first + ".ini";
    process.start((command + " -write_ini " + tmp_file).toQString());
    process.waitForFinished();
    Param p;
    ParamXMLFile pf;
    pf.load(tmp_file, p);
    File::remove(tmp_file);
    ofstream f_html((String("output/") + prefix + it->first + ".html").c_str());
    convertINI2HTML(p, f_html);
    f_html.close();
  }
  return errors_occured;
}

int main(int argc, char ** argv)
{
  if (argc != 2)
  {
    cerr << "Please specify the path where the TOPP/UTIL binaries are located." << endl;
    return EXIT_FAILURE;
  }
  
  String binary_directory = String(argv[1]).ensureLastChar('/');

  if(!File::exists(binary_directory))
  {
    cerr << "The given binary directory does not exist. Aborting." << endl;
    return EXIT_FAILURE;
  }
  
  //TOPP tools
  ToolListType topp_tools = ToolHandler::getTOPPToolList(true);   // include GenericWrapper (can be called with --help without error, even though it has a type)
  topp_tools["TOPPView"] = Internal::ToolDescription();   // these two need to be excluded from writing an INI file later!
  topp_tools["TOPPAS"] = Internal::ToolDescription();
  //UTILS
  ToolListType util_tools = ToolHandler::getUtilList();

  bool errors_occured = generate(topp_tools, "TOPP_", binary_directory) || generate(util_tools, "UTILS_", binary_directory);

  if (errors_occured)
  {
    // errors occurred while generating the TOPP CLI docu .. tell the user
    cerr << "Errors occurred while generating the command line documentation for some of the TOPP tools/UTILS." << endl;
    return EXIT_FAILURE;
  }
  else
  {
    return EXIT_SUCCESS;
  }
}
