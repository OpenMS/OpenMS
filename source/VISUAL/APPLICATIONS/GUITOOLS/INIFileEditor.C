// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/VISUAL/APPLICATIONS/INIFileEditorWindow.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

// Qt
#include <QtGui/QStyleFactory>
#include <QTextCodec>


#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif


using namespace OpenMS;
using namespace std;

/**
    @page TOPP_INIFileEditor INIFileEditor

    @brief Can be used to visually edit INI files of TOPP tools.

    The values can be edited by double-clicking or pressing F2.

    The documentation of each value is shown in the text area on the bottom of the widget.

    @image html INIFileEditor.png
*/

int main(int argc, const char** argv)
{
#if  defined(__APPLE__)
  // we do not want to load plugins as this leads to serious problems
  // when shipping on mac os x
  QApplication::setLibraryPaths(QStringList());
#endif
  
  // ensure correct encoding of paths
  QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
  
  Map<String, String> option_lists;
  Map<String, String> options;
  options["-print"] = "print";
  Map<String, String> flags;
  flags["--help"] = "help";
  Param param;
  param.parseCommandLine(argc, argv, options, flags, option_lists);

  //catch command line errors
  if (param.exists("help") //help requested
     || argc > 3 //too many arguments
     || (argc == 3 && !param.exists("print")) //three argument but no -print
     || (param.exists("print") && param.getValue("print") == "") //-print but no file given
      )
  {
    cerr << endl
         << "INIFileEditor -- An editor for OpenMS configuration files." << endl
         << endl
         << "Usage:" << endl
         << " INIFileEditor [options] [file]" << endl
         << endl
         << "Options are:" << endl
         << " --help         Shows this help and exits" << endl
         << " -print <file>  Prints the content of the file to the command line and exits" << endl
         << endl;
    return 0;
  }

  //print a ini file as text
  if (param.exists("print"))
  {
    Param data;
    ParamXMLFile paramFile;
    try
    {
      paramFile.load(param.getValue("print"), data);
      for (Param::ParamIterator it = data.begin(); it != data.end(); ++it)
      {
        cout << it.getName() << " = " << it->value << endl;
      }
    }
    catch (Exception::BaseException& e)
    {
      LOG_ERROR << "Error while parsing file '" << param.getValue("print") << "'\n";
      LOG_ERROR << e << "\n";
    }

    return 0;
  }

  //Create window
  QApplicationTOPP app(argc, const_cast<char**>(argv));

  //set plastique style unless windows / mac style is available
  if (QStyleFactory::keys().contains("windowsxp", Qt::CaseInsensitive))
  {
    app.setStyle("windowsxp");
  }
  else if (QStyleFactory::keys().contains("macintosh", Qt::CaseInsensitive))
  {
    app.setStyle("macintosh");
  }
  else if (QStyleFactory::keys().contains("plastique", Qt::CaseInsensitive))
  {
    app.setStyle("plastique");
  }

  INIFileEditorWindow editor_window;

  //Open passed file
  if (argc == 2)
  {
    //cout << "OPEN: "  << argv[1] << endl;
    editor_window.openFile(argv[1]);
  }

#ifdef OPENMS_WINDOWSPLATFORM
  FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
  AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

  editor_window.show();
  return app.exec();
}
