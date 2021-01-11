// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

/**
  @page TOPP_TOPPView TOPPView

  TOPPView is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzML, mzData, mzXML
  and several other file formats. It also supports viewing data from an %OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC-MS map and a MS raw spectrum:

  @image html TOPPView.png

  More information about TOPPView can be found in the @ref TOPP_tutorial.

  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TOPPView.cli
*/

//QT
#include <QtWidgets/QSplashScreen>
#include <QMessageBox>

//OpenMS
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/SYSTEM/StopWatch.h>


using namespace OpenMS;
using namespace std;

//STL
#include <iostream>
#include <map>
#include <vector>

#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "TOPPView";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
  cerr << endl
       << tool_name << " -- A viewer for mass spectrometry data." << endl
       << endl
       << "Usage:" << endl
       << " " << tool_name << " [options] [files]" << endl
       << endl
       << "Options are:" << endl
       << "  --help           Shows this help" << endl
       << "  -ini <File>      Sets the INI file (default: ~/.TOPPView.ini)" << endl
       << endl
       << "Hints:" << endl
       << " - To open several files in one window put a '+' in between the files." << endl
       << " - '@bw' after a map file displays the dots in a white to black gradient." << endl
       << " - '@bg' after a map file displays the dots in a grey to black gradient." << endl
       << " - '@b'  after a map file displays the dots in black." << endl
       << " - '@r'  after a map file displays the dots in red." << endl
       << " - '@g'  after a map file displays the dots in green." << endl
       << " - '@m'  after a map file displays the dots in magenta." << endl
       << " - Example: '" << tool_name << " 1.mzML + 2.mzML @bw + 3.mzML @bg'" << endl
       << endl;
}

int main(int argc, const char** argv)
{
  //list of all the valid options
  Map<String, String> valid_options, valid_flags, option_lists;
  valid_flags["--help"] = "help";
  valid_options["-ini"] = "ini";

  Param param;
  param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

  // '--help' given
  if (param.exists("help"))
  {
    print_usage();
    return 0;
  }

  // test if unknown options were given
  if (param.exists("unknown"))
  {
    // if TOPPView is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if (!(param.getValue("unknown").toString().hasSubstring("-psn") && !param.getValue("unknown").toString().hasSubstring(", ")))
    {
      cout << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage();
      return 1;
    }
  }

  try
  {
    QApplicationTOPP a(argc, const_cast<char**>(argv));
    a.connect(&a, &QApplicationTOPP::lastWindowClosed, &a, &QApplicationTOPP::quit);

    TOPPViewBase tb;
    a.connect(&a, &QApplicationTOPP::fileOpen, &tb, &TOPPViewBase::openFile);
    tb.show();

    // Create the splashscreen that is displayed while the application loads (version is drawn dynamically)
    QPixmap qpm(":/TOPPView_Splashscreen.png");
    QPainter pt_ver(&qpm);
    pt_ver.setFont(QFont("Helvetica [Cronyx]", 15, 2, true));
    pt_ver.setPen(QColor(44, 50, 152));
    pt_ver.drawText(490, 94, VersionInfo::getVersion().toQString());
    QSplashScreen splash_screen(qpm);
    splash_screen.show();

    QApplication::processEvents();
    StopWatch stop_watch;
    stop_watch.start();

    if (param.exists("ini"))
    {
      tb.loadPreferences((String)param.getValue("ini"));
    }

    //load command line files
    if (param.exists("misc"))
    {
      tb.loadFiles(param.getValue("misc"), &splash_screen);
    }

    // We are about to show the application.
    // Proper time to remove the splashscreen, if at least 1.5 seconds have passed...
    while (stop_watch.getClockTime() < 1.5) /*wait*/
    {
    }
    stop_watch.stop();
    splash_screen.close();

#ifdef OPENMS_WINDOWSPLATFORM
    FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
    AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif
    return a.exec();
  }
  //######################## ERROR HANDLING #################################
  catch (Exception::UnableToCreateFile& e)
  {
    cout << String("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotFound& e)
  {
    cout << String("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileNotReadable& e)
  {
    cout << String("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::FileEmpty& e)
  {
    cout << String("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::ParseError& e)
  {
    cout << String("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::InvalidValue& e)
  {
    cout << String("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }
  catch (Exception::BaseException& e)
  {
    cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
  }

  return 1;
}
