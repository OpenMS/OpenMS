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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/VISUAL/APPLICATIONS/IDEvaluationBase.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

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
    @page UTILS_IDEvaluation IDEvaluation

    @brief Computes a 'q-value vs. #PSM' plot to visualize the number identifications for a certain q-value.


  An arbitrary number of idXML files resulting from a target+decoy search can be provided as input.

  Since the q-value can be computed independently from a scoring scheme, no further preprocessing (like IDPep or FDR)
  is required, apart from a target-decoy annotation! I.e., apply PeptideIndexer to the immediate output of a search engine
  (or ConsensusID) and use this as input to this tool.


    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_IDEvaluation.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_IDEvaluation.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDEvaluation :
  public TOPPBase
{
public:
  TOPPIDEvaluation() :
    TOPPBase("IDEvaluation",
             "Computes a 'q-value vs. #PSM' plot to visualize the number identifications for a certain q-value.", false)
  {
    char * dummy = "dummy"; int argc = 1;  QApplication a(argc, &dummy);
    out_formats_ = IDEvaluationBase().getSupportedImageFormats(); // can only be called if a QApplication is present...
  }

protected:
  StringList out_formats_; //< valid output formats for image

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param p_my;

    Param p = FalseDiscoveryRate().getDefaults();
    p_my.insert("fdr:", p.copy("use_all_hits"));

    char * dummy = "dummy"; int argc = 1;  QApplication a(argc, &dummy);
    p_my.insert("image:", IDEvaluationBase().getParameters().copy("image:", true)); // can only be called if a QApplication is present...
    return p_my;
  }

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<file>", StringList::create(""), "Input file(s)", false);
    setValidFormats_("in", StringList::create("idXML"));


    registerOutputFile_("out", "<file>", "", "Output file (if given, no GUI will be displayed)", false);
    setValidFormats_("out", out_formats_, false);
    registerStringOption_("out_type", "<file type>", "", "The image format. Set this if you want to force a format not reflected by the 'out' filename.", false);
    setValidStrings_("out_type", out_formats_);
    registerOutputFile_("out_csv", "<file>", "", "Optional output of points as table for manual post-processing.", false);

    registerDoubleOption_("q_min", "<float>", 0.0, "Minimal q-value in plot.", false);
    setMinFloat_("q_min", 0.0);
    setMaxFloat_("q_min", 1.0);
    registerDoubleOption_("q_max", "<float>", 0.4, "Maximal q-value in plot.", false);
    setMinFloat_("q_max", 0.0);
    setMaxFloat_("q_max", 1.0);

    registerSubsection_("algorithm", "Additional parameters for FDR and image sizes.");
  }

  ExitCodes main_(int argc, const char ** argv)
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    StringList in_list = getStringList_("in");
    String out = getStringOption_("out").trim();
    String format = getStringOption_("out_type").trim();
    if (out != "" && format == "") // get from filename
    {
      try
      {
        format = out.suffix('.');
      }
      catch (Exception::ElementNotFound & /*e*/)
      {
        format = "nosuffix";
      }
      // check if format is valid:
      if (!out_formats_.contains(format.toLower()))
      {
        LOG_ERROR << "No explicit image output format was provided via 'out_type', and the suffix ('" << format << "') does not resemble a valid type. Please fix one of them." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    DoubleReal q_min = getDoubleOption_("q_min");
    DoubleReal q_max = getDoubleOption_("q_max");
    if (q_min >= q_max)
    {
      LOG_ERROR << "The parameter 'q_min' must be smaller than 'q_max'. Quitting..." << std::endl;
      return ILLEGAL_PARAMETERS;
    }


    QApplication a(argc, const_cast<char **>(argv));

    IDEvaluationBase * mw = new IDEvaluationBase();
    Param alg_param = mw->getParameters();
    alg_param.insert("", getParam_().copy("algorithm:", true));
    mw->setParameters(alg_param);
    mw->loadFiles(in_list);
    mw->setVisibleArea(q_min, q_max);
    mw->show(); // required to get the size of the images right

    if (out != "") // save as image and exit
    {
      String error;
      bool r = mw->exportAsImage(out.toQString(), error, format.toQString());
      if (r) return EXECUTION_OK;
      else
      {
        LOG_ERROR << error << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }

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
  TOPPIDEvaluation tool;
  return tool.main(argc, argv);
}

/// @endcond
