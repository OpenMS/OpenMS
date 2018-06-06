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
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/MATH/MISC/BilinearInterpolation.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/MultiGradient.h>
#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>


#include <QtGui/QImage>
#include <QtGui/QPainter>


using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_ImageCreator ImageCreator

    @brief Transforms an LC-MS map into a png image.

    The input is first resampled into a matrix using
    bilinear forward resampling.  Then the content of the matrix is written to
    an image file.  The output has a uniform spacing in both dimensions regardless
    of the input.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_ImageCreator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_ImageCreator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPImageCreator :
  public TOPPBase
{
public:
  TOPPImageCreator() :
    TOPPBase("ImageCreator",
             "Transforms an LC-MS map into an image.", false)
  {
    out_formats_ = ListUtils::create<String>("png,jpg,bmp,tiff,ppm");
  }

protected:
  StringList out_formats_; ///< valid output formats for image

  void addPoint_(int x, int y, QImage & image, QColor color = Qt::black,
                 Size size = 2)
  {
    int h = image.height(), w = image.width();
    vector<int> xs(1, x), ys(1, y);
    if (size == 2)
    {
      int xtemp[] = {x - 1, x, x, x + 1};
      int ytemp[] = {y, y - 1, y + 1, y};
      xs = vector<int>(xtemp, xtemp + 4);
      ys = vector<int>(ytemp, ytemp + 4);
    }
    else if (size == 3)
    {
      int xtemp[] = {x - 2, x - 1, x - 1, x, x, x + 1, x + 1, x + 2};
      int ytemp[] = {y, y + 1, y - 1, y + 2, y - 2, y + 1, y - 1, y};
      xs = vector<int>(xtemp, xtemp + 8);
      ys = vector<int>(ytemp, ytemp + 8);
    }
    for (Size i = 0; i < xs.size(); ++i)
    {
      int xi = xs[i], yi = ys[i];
      if ((xi > 0) && (xi < w) && (yi > 0) && (yi < h))
      {
        image.setPixel(xi, yi, color.rgb());
      }
    }
  }

  void addFeatureBox_(int lower_mz, int lower_rt, int upper_mz, int upper_rt, QImage & image, QColor color = Qt::black)
  {
    QPainter * painter = new QPainter(&image);
    painter->setPen(color);
    painter->drawRect(QRect(lower_rt, lower_mz, upper_rt - lower_rt, upper_mz - lower_mz));
    delete painter;
  }

  void markMS2Locations_(PeakMap & exp, QImage & image, bool transpose,
                         QColor color, Size size)
  {
    double xcoef = image.width(), ycoef = image.height();
    if (transpose)
    {
      xcoef /= exp.getMaxRT() - exp.getMinRT();
      ycoef /= exp.getMaxMZ() - exp.getMinMZ();
    }
    else
    {
      xcoef /= exp.getMaxMZ() - exp.getMinMZ();
      ycoef /= exp.getMaxRT() - exp.getMinRT();
    }
    for (PeakMap::Iterator spec_iter = exp.begin();
         spec_iter != exp.end(); ++spec_iter)
    {
      if (spec_iter->getMSLevel() == 2)
      {
        double mz = spec_iter->getPrecursors()[0].getMZ();
        double rt = exp.getPrecursorSpectrum(spec_iter)->getRT();
        int x, y;
        if (transpose)
        {
          x = int(xcoef * (rt - exp.getMinRT()));
          y = int(ycoef * (exp.getMaxMZ() - mz));
        }
        else
        {
          x = int(xcoef * (mz - exp.getMinMZ()));
          y = int(ycoef * (exp.getMaxRT() - rt));
        }
        addPoint_(x, y, image, color, size);  //mark MS2
      }
    }
  }

  void markFeatureLocations_(FeatureMap & feature_map, PeakMap & exp, QImage & image, bool transpose, QColor color)
  {
    double xcoef = image.width(), ycoef = image.height();
    if (transpose)
    {
      xcoef /= exp.getMaxRT() - exp.getMinRT();
      ycoef /= exp.getMaxMZ() - exp.getMinMZ();
    }
    else
    {
      xcoef /= exp.getMaxMZ() - exp.getMinMZ();
      ycoef /= exp.getMaxRT() - exp.getMinRT();
    }

    for (FeatureMap::Iterator feat_iter = feature_map.begin();
         feat_iter != feature_map.end(); ++feat_iter)
    {
      const ConvexHull2D convex_hull = feat_iter->getConvexHull();
      DBoundingBox<2> box = convex_hull.getBoundingBox();
      double rt = feat_iter->getRT();
      double mz = feat_iter->getMZ();
      double lower_mz = box.minY();
      double lower_rt = box.minX();
      double upper_mz = box.maxY();
      double upper_rt = box.maxX();

      int lx, ly, ux, uy, cx, cy;
      if (transpose)
      {
        lx = int(xcoef * (lower_rt - exp.getMinRT()));
        ly = int(ycoef * (exp.getMaxMZ() - lower_mz));
        ux = int(xcoef * (upper_rt - exp.getMinRT()));
        uy = int(ycoef * (exp.getMaxMZ() - upper_mz));
        cx = int(xcoef * (rt - exp.getMinRT()));
        cy = int(ycoef * (mz - lower_mz));
      }
      else
      {
        lx = int(xcoef * (lower_mz - exp.getMinMZ()));
        ly = int(ycoef * (exp.getMaxRT() - lower_rt));
        ux = int(xcoef * (upper_mz - exp.getMinMZ()));
        uy = int(ycoef * (exp.getMaxRT() - upper_rt));
        cx = int(xcoef * (mz - exp.getMinMZ()));
        cy = int(ycoef * (exp.getMaxRT() - rt));
      }

      addFeatureBox_(ly, lx, uy, ux, image, color);
      addPoint_(cx, cy, image, Qt::black); // mark center
    }
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("in_featureXML", "<file>", "", "input file ", false);
    setValidFormats_("in_featureXML", ListUtils::create<String>("featureXML"));

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", out_formats_, false);
    registerStringOption_("out_type", "<file type>", "", "The image format. Set this if you want to force a format not reflected by the 'out' filename.", false);
    setValidStrings_("out_type", out_formats_);

    registerIntOption_("width", "<number>", 1024, "Number of pixels in m/z dimension.\nIf 0, one pixel per Th.", false);
    setMinInt_("width", 0);
    registerIntOption_("height", "<number>", 1024, "Number of pixels in RT dimension.\nIf 0, one pixel per spectrum.", false);
    setMinInt_("height", 0);
    registerStringOption_("background_color", "<color>", "#FFFFFF", "Background color e.g.: \"#FF0000\" to choose red as background color", false);
    registerStringOption_("feature_color", "<color>", "#000000", "Feature color e.g.: \"#00FF00\" to choose green as feature color", false);

    registerStringOption_("gradient", "<gradient>", "", "Intensity gradient that defines colors for the range between 0 and 100.\n"
                                                        "Example: '0,#FFFFFF;50,#FF0000;100,#000000'", false);
    registerDoubleOption_("max_intensity", "<int>", 0, "Maximum peak intensity used to determine range for colors.\n"
                                                       "If 0, this is determined from the data.", false);
    registerFlag_("log_intensity", "Apply logarithm to intensity values");
    registerFlag_("transpose", "flag to transpose the resampled matrix (RT vs. m/z).\n"
                               "Per default, dimensions run bottom-up in RT and left-right in m/z.");
    registerFlag_("precursors", "Mark locations of MS2 precursors.\n");
    registerStringOption_("precursor_color", "<color>", "#000000", "Color for precursor marks (color code or word, e.g. 'black') (requires 'precursors' flag to be active)", false);
    registerIntOption_("precursor_size", "<number>", 2,
                       "Size of the precursor marks (requires 'precursors' flag to be active)", false);
    setMinInt_("precursor_size", 1);
    setMaxInt_("precursor_size", 3);
  }

  ExitCodes main_(int, const char **) override
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    String in = getStringOption_("in");
    String in_featureXML = getStringOption_("in_featureXML");
    String out = getStringOption_("out");
    String format = getStringOption_("out_type");
    if (format.trim() == "") // get from filename
    {
      try
      {
        format = out.suffix('.');
      }
      catch (Exception::ElementNotFound & /*e*/)
      {
        format = "nosuffix";
      }
      StringListUtils::toUpper(out_formats_);
      if (!ListUtils::contains(out_formats_, format.toUpper()))
      {
        LOG_ERROR << "No explicit image output format was provided via 'out_type', and the suffix ('" << format << "') does not resemble a valid type. Please fix one of them." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }
    PeakMap exp;
    MzMLFile f;
    f.setLogType(log_type_);
    f.load(in, exp);

    exp.updateRanges(1);

    SignedSize rows = getIntOption_("height");
    if (rows == 0)
    {
      rows = exp.size();
    }
    if (rows <= 0)
    {
      writeLog_("Error: Zero rows is not possible.");
      return ILLEGAL_PARAMETERS;
    }

    SignedSize cols = getIntOption_("width");
    if (cols == 0)
    {
      cols = UInt(ceil(exp.getMaxMZ() - exp.getMinMZ()));
    }
    if (cols <= 0)
    {
      writeLog_("Error: Zero columns is not possible.");
      return ILLEGAL_PARAMETERS;
    }

    //----------------------------------------------------------------
    //Do the actual resampling
    BilinearInterpolation<double, double> bilip;
    bilip.getData().resize(rows, cols);

    if (!getFlag_("transpose"))
    {
      // scans run bottom-up:
      bilip.setMapping_0(0, exp.getMaxRT(), rows - 1, exp.getMinRT());
      // peaks run left-right:
      bilip.setMapping_1(0, exp.getMinMZ(), cols - 1, exp.getMaxMZ());

      for (PeakMap::Iterator spec_iter = exp.begin();
           spec_iter != exp.end(); ++spec_iter)
      {
        if (spec_iter->getMSLevel() != 1) continue;
        for (PeakMap::SpectrumType::ConstIterator peak1_iter =
               spec_iter->begin(); peak1_iter != spec_iter->end();
             ++peak1_iter)
        {
          bilip.addValue(spec_iter->getRT(), peak1_iter->getMZ(),
                         peak1_iter->getIntensity());
        }
      }
    }
    else     // transpose
    {
      // spectra run bottom-up:
      bilip.setMapping_0(0, exp.getMaxMZ(), rows - 1, exp.getMinMZ());
      // scans run left-right:
      bilip.setMapping_1(0, exp.getMinRT(), cols - 1, exp.getMaxRT());

      for (PeakMap::Iterator spec_iter = exp.begin();
           spec_iter != exp.end(); ++spec_iter)
      {
        if (spec_iter->getMSLevel() != 1) continue;
        for (PeakMap::SpectrumType::ConstIterator peak1_iter =
               spec_iter->begin(); peak1_iter != spec_iter->end();
             ++peak1_iter)
        {
          bilip.addValue(peak1_iter->getMZ(), spec_iter->getRT(),
                         peak1_iter->getIntensity());
        }
      }
    }

    //----------------------------------------------------------------
    //create and store image
    int scans = (int) bilip.getData().sizePair().first;
    int peaks = (int) bilip.getData().sizePair().second;

    MultiGradient gradient;
    String gradient_str = getStringOption_("gradient");
    if (gradient_str != "")
    {
      gradient.fromString(String("Linear|") + gradient_str);
    }
    else
    {
      gradient.fromString("Linear|0,#FFFFFF;2,#FFFF00;11,#FFAA00;32,#FF0000;55,#AA00FF;78,#5500FF;100,#000000");
    }

    bool use_log = getFlag_("log_intensity");
    writeDebug_("log_intensity: " + String(use_log), 1);

    QImage image(peaks, scans, QImage::Format_RGB32);
    string s = getStringOption_("background_color");
    QColor background_color(s.c_str());

    string feature_color_string = getStringOption_("feature_color");
    QColor feature_color(feature_color_string.c_str());

    QPainter * painter = new QPainter(&image);
    painter->setPen(background_color);
    painter->fillRect(0, 0, peaks, scans, Qt::SolidPattern);
    delete painter;

    double factor = getDoubleOption_("max_intensity");
    if (factor == 0)
    {
      factor = (*std::max_element(bilip.getData().begin(), bilip.getData().end()));
    }
    // logarithmize max. intensity as well:
    if (use_log) factor = std::log(factor);

    factor /= 100.0;
    for (int i = 0; i < scans; ++i)
    {
      for (int j = 0; j < peaks; ++j)
      {
        double value = bilip.getData().getValue(i, j);
        if (use_log) value = std::log(value);
        if (value > 1e-4)
        {
          image.setPixel(j, i, gradient.interpolatedColorAt(value / factor).rgb());
        }
        else
        {
          image.setPixel(j, i, background_color.rgb());
        }
      }
    }

    if (getFlag_("precursors"))
    {
      markMS2Locations_(exp, image, getFlag_("transpose"),
                        getStringOption_("precursor_color").toQString(),
                        Size(getIntOption_("precursor_size")));
    }

    if (!in_featureXML.empty())
    {
      FeatureMap feature_map;
      FeatureXMLFile ff;
      ff.load(in_featureXML, feature_map);
      markFeatureLocations_(feature_map, exp, image, getFlag_("transpose"), feature_color);
    }

    if (image.save(out.toQString(), format.c_str())) return EXECUTION_OK;
    else return CANNOT_WRITE_OUTPUT_FILE;
  }

};


int main(int argc, const char ** argv)
{
  TOPPImageCreator tool;
  return tool.main(argc, argv);
}

/// @endcond
