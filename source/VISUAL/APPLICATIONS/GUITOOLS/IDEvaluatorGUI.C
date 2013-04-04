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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/IDEvaluationBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>

#include <QtGui/QImage>
#include <QPainter>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QApplication>

#ifdef OPENMS_WINDOWSPLATFORM
#   ifndef _WIN32_WINNT
#       define _WIN32_WINNT 0x0501 // Win XP (and above)
#   endif
#   include <Windows.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_IDEvaluatorGUI IDEvaluatorGUI

    @brief Computes a 'q-value vs. #PSM' plot to visualize the number identifications for a certain q-value.


  An arbitrary number of idXML files resulting from a target+decoy search can be provided as input.

  Since the q-value can be computed independently from a scoring scheme, no further preprocessing (like IDPep or FDR)
  is required, apart from a target-decoy annotation! I.e., apply PeptideIndexer to the immediate output of a search engine
  (or ConsensusID) and use this as input to this tool.


  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_IDEvaluatorGUI.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_IDEvaluatorGUI.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDEvaluatorGUI :
  public TOPPBase
{
public:
  TOPPIDEvaluatorGUI() :
    TOPPBase("IDEvaluatorGUI",
             "Computes a 'q-value vs. #PSM' plot to visualize the number identifications for a certain q-value.", false, false, false)
  {
    int argc = 1;
    const char* c = "IDEvaluatorGUI";
    const char** argv = &c;
    QApplication a(argc, const_cast<char**>(argv));
    out_formats_ = IDEvaluationBase().getSupportedImageFormats(); // can only be called if a QApplication is present...
  }

protected:
  StringList out_formats_; //< valid output formats for image

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param p_my;

    Param p = FalseDiscoveryRate().getDefaults();
    p_my.insert("fdr:", p.copy("use_all_hits"));

    int argc = 1;
    const char* c = "IDEvaluatorGUI";
    const char** argv = &c;
    QApplication a(argc, const_cast<char**>(argv));
    p_my.insert("image:", IDEvaluationBase().getParameters().copy("image:", true)); // can only be called if a QApplication is present...
    return p_my;
  }

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<file>", StringList::create(""), "Input file(s)", false);
    setValidFormats_("in", StringList::create("idXML"));

    registerSubsection_("algorithm", "Additional parameters for FDR and image sizes.");
  }

  ExitCodes main_(int argc, const char ** argv)
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    StringList in_list = getStringList_("in");

    QApplicationTOPP a(argc, const_cast<char **>(argv));

    IDEvaluationBase * mw = new IDEvaluationBase();
    Param alg_param = mw->getParameters();
    alg_param.insert("", getParam_().copy("algorithm:", true));
    mw->setParameters(alg_param);
    if (!mw->loadFiles(in_list))
    {
      LOG_ERROR << "Tool failed. See above." << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    };

    mw->show();

#ifdef OPENMS_WINDOWSPLATFORM
    FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
    AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

    int result = a.exec();
    delete(mw);
    if (result) return UNKNOWN_ERROR;
    else return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPIDEvaluatorGUI tool;
  return tool.main(argc, argv);
}

/// @endcond
