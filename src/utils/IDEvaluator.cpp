// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/VISUAL/APPLICATIONS/IDEvaluationBase.h>
#include <OpenMS/VISUAL/APPLICATIONS/MISC/QApplicationTOPP.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

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
  @page UTILS_IDEvaluator IDEvaluator

  @brief Computes a 'q-value vs. #PSM' plot which is saved as an image to visualize the number identifications for a certain q-value.


  An arbitrary number of idXML files resulting from a target+decoy search can be provided as input.

  Since the q-value can be computed independently from a scoring scheme, no further preprocessing (like IDPep or FDR)
  is required, apart from a target-decoy annotation! I.e., apply PeptideIndexer to the immediate output of a search engine
  (or ConsensusID) and use this as input to this tool.


  @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_IDEvaluator.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_IDEvaluator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDEvaluator :
  public TOPPBase
{
public:
  TOPPIDEvaluator() :
    TOPPBase("IDEvaluator",
             "Computes a 'q-value vs. #PSM' plot which is saved as an image to visualize the number identifications for a certain q-value.", false)
  {
    out_formats_ = IDEvaluationBase().getSupportedImageFormats(); // can only be called if a QApplication is present...
  }

protected:
  StringList out_formats_; ///< valid output formats for image

  Param getSubsectionDefaults_(const String& /*section*/) const override
  {
    Param p_my;

    Param p = FalseDiscoveryRate().getDefaults();
    p_my.insert("fdr:", p.copy("use_all_hits"));

    p_my.insert("image:", IDEvaluationBase().getParameters().copy("image:", true)); // can only be called if a QApplication is present...
    return p_my;
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file>", ListUtils::create<String>(""), "Input file(s)", false);
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file (if given, no GUI will be displayed)", false);
    setValidFormats_("out", out_formats_, false);
    registerStringOption_("out_type", "<file type>", "", "The image format. Set this if you want to force a format not reflected by the 'out' filename.", false);
    setValidStrings_("out_type", out_formats_);
    registerOutputFile_("out_csv", "<file>", "", "Optional output of points as table for manual post-processing.", false);
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));

    registerDoubleOption_("q_min", "<float>", 0.0, "Minimal q-value in plot.", false);
    setMinFloat_("q_min", 0.0);
    setMaxFloat_("q_min", 1.0);
    registerDoubleOption_("q_max", "<float>", 0.4, "Maximal q-value in plot.", false);
    setMinFloat_("q_max", 0.0);
    setMaxFloat_("q_max", 1.0);

    registerSubsection_("algorithm", "Additional parameters for FDR and image sizes.");
  }

  ExitCodes main_(int, const char**) override
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    StringList in_list = getStringList_("in");
    String out = getStringOption_("out");
    String out_csv = getStringOption_("out_csv");
    String format = getStringOption_("out_type");

    if (out.empty() && out_csv.empty())
    {
      LOG_ERROR << "Neither 'out' nor 'out_csv' were provided. Please assign at least one of them." << std::endl;
      return ILLEGAL_PARAMETERS;
    }

    if (!out.empty() && format == "") // get from filename
    {
      try
      {
        format = out.suffix('.');
      }
      catch (Exception::ElementNotFound& /*e*/)
      {
        format = "nosuffix";
      }
      // check if format is valid:
      if (!ListUtils::contains(out_formats_, format.toLower()))
      {
        LOG_ERROR << "No explicit image output format was provided via 'out_type', and the suffix ('" << format << "') does not resemble a valid type. Please fix one of them." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }

    double q_min = getDoubleOption_("q_min");
    double q_max = getDoubleOption_("q_max");
    if (q_min >= q_max)
    {
      LOG_ERROR << "The parameter 'q_min' must be smaller than 'q_max'. Quitting..." << std::endl;
      return ILLEGAL_PARAMETERS;
    }


    IDEvaluationBase* mw = new IDEvaluationBase();
    Param alg_param = mw->getParameters();
    alg_param.insert("", getParam_().copy("algorithm:", true));
    mw->setParameters(alg_param);

    if (!mw->loadFiles(in_list))
    {
      LOG_ERROR << "Tool failed. See above." << std::endl;
      return INCOMPATIBLE_INPUT_DATA;
    }
    mw->setVisibleArea(q_min, q_max);

    if (!out.empty()) // save as image and exit
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

    if (!out_csv.empty())
    {
      TextFile tf;
      for (Size i = 0; i < mw->getPoints().size(); ++i)
      {
        MSSpectrum s = mw->getPoints()[i];
        StringList sl1;
        StringList sl2;
        for (Size j = 0; j < s.size(); ++j)
        {
          sl1.push_back(s[j].getMZ());
          sl2.push_back(s[j].getIntensity());
        }
        tf.addLine(String("# ") + String(s.getMetaValue("search_engine")));
        tf.addLine(ListUtils::concatenate(sl1, ","));
        tf.addLine(ListUtils::concatenate(sl2, ","));
      }
      tf.store(out_csv);
    }

    delete(mw);
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{

  /*
   * Initializing Qt Application which gets its own version of argc and argv
   * (apparently, it will modify argv and does not use a const version of it).
   *
   * We need to make sure to only initialize QApplication since it is a
   * singleton. Creating multiple ones will cause a segfault with Qt 4.8.1
   * and possibly other versions as well.
   *
   * See bug 569
   *
  */
  int argc_ = 1;
  const char* c = "IDEvaluator";
  const char** argv_ = &c;
  QApplication app(argc_, const_cast<char**>(argv_)); // no QApplicationTOPP, since any exception thrown will result in an abort of a cmd line tool anyways (no real GUI here)

  TOPPIDEvaluator tool;
  return tool.main(argc, argv);
}

/// @endcond
