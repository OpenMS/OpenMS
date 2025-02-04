// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/PROCESSING/RESAMPLING/LinearResampler.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/ML/INTERPOLATION/BilinearInterpolation.h>
#include <OpenMS/VISUAL/MultiGradient.h>


#include <QtGui/QImage>
#include <QtGui/QPainter>


using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ImageCreator ImageCreator

@brief Transforms an LC-MS map into a png image.

The input is first resampled into a matrix using bilinear forward resampling.
Then the content of the matrix is written to an image file.
The output has a uniform spacing in both dimensions regardless of the input.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ImageCreator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ImageCreator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPImageCreator :
  public TOPPBase
{
public:
  TOPPImageCreator() :
    TOPPBase("ImageCreator",
             "Transforms an LC-MS map into an image.", false), 
    out_formats_({"png", "jpg", "bmp", "tiff", "ppm"}) // all in lower case!
  {
  }

protected:
  const StringList out_formats_; ///< valid output formats for image

  void addPoint_(int x, int y, QImage& image, QColor color = Qt::black,
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

  void addFeatureBox_(int lower_mz, int lower_rt, int upper_mz, int upper_rt, QImage& image, QColor color = Qt::black)
  {
    QPainter * painter = new QPainter(&image);
    painter->setPen(color);
    painter->drawRect(QRect(lower_rt, lower_mz, upper_rt - lower_rt, upper_mz - lower_mz));
    delete painter;
  }

  void markMS2Locations_(PeakMap& exp, QImage& image, bool transpose,
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

  void markFeatureLocations_(FeatureMap& feature_map, PeakMap& exp, QImage& image, bool transpose, QColor color)
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
    setValidFormats_("in", {"mzML"});
    registerInputFile_("in_featureXML", "<file>", "", "input file ", false);
    setValidFormats_("in_featureXML", {"featureXML"});

    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", out_formats_, false);
    registerStringOption_("out_type", "<file type>", "", "The image format. Set this if you want to force a format not reflected by the 'out' filename.", false);
    setValidStrings_("out_type", out_formats_);

    registerStringOption_("rt", "[min]:[max]", ":", "Retention time range to extract", false);
    registerStringOption_("mz", "[min]:[max]", ":", "Mass-to-charge range to extract", false);

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
    registerFlag_("transpose", "Flag to transpose the resampled matrix (RT vs. m/z).\n"
                               "Per default, dimensions run bottom-up in RT and left-right in m/z.");
    registerFlag_("precursors", "Mark locations of MS2 precursors.\n");
    registerStringOption_("precursor_color", "<color>", "#000000", "Color for precursor marks (color code or word, e.g. 'black') (requires 'precursors' flag to be active)", false);
    registerIntOption_("precursor_size", "<number>", 2,
                       "Size of the precursor marks (requires 'precursors' flag to be active)", false);
    setMinInt_("precursor_size", 1);
    setMaxInt_("precursor_size", 3);
  }

  ExitCodes main_(int, const char**) override
  {
    //----------------------------------------------------------------
    // load data
    //----------------------------------------------------------------
    String in = getStringOption_("in");
    String in_featureXML = getStringOption_("in_featureXML");
    String out = getStringOption_("out");
    String format = getStringOption_("out_type");
    if (format.trim().empty()) // get from filename
    {
      try
      {
        format = out.suffix('.');
      }
      catch (Exception::ElementNotFound& /*e*/)
      {
        format = "nosuffix";
      }
      if (!ListUtils::contains(out_formats_, format.toLower()))
      {
        OPENMS_LOG_ERROR << "No explicit image output format was provided via 'out_type', and the suffix ('" << format << "') does not resemble a valid type. Please fix one of them." << std::endl;
        return ILLEGAL_PARAMETERS;
      }
    }
    const double init = numeric_limits<double>::max();
    double rt_min = -init, rt_max = init, mz_min = -init, mz_max = init;
    bool filter_rt = parseRange_(getStringOption_("rt"), rt_min, rt_max);
    if (rt_min > rt_max) swap(rt_min, rt_max);
    bool filter_mz = parseRange_(getStringOption_("mz"), mz_min, mz_max);
    if (mz_min > mz_max) swap(mz_min, mz_max);
    bool show_precursors = getFlag_("precursors");

    PeakMap exp;
    FileHandler f;
    if (filter_rt) f.getOptions().setRTRange(DRange<1>(rt_min, rt_max));
    if (filter_mz) f.getOptions().setMZRange(DRange<1>(mz_min, mz_max));
    if (!show_precursors) f.getOptions().setMSLevels({1});
    f.loadExperiment(in, exp, {FileTypes::MZML}, log_type_);
    if (filter_mz && show_precursors)
    {
      // MS2 spectra were not filtered by precursor m/z, remove them now:
      auto predicate =
        InPrecursorMZRange<MSSpectrum>(mz_min, mz_max, true);
      exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), predicate),
                             exp.end());
    }
    exp.updateRanges(1);

    Size rows = getIntOption_("height"), cols = getIntOption_("width");
    if (rows == 0) rows = exp.size();
    if (cols == 0) cols = UInt(ceil(exp.getMaxMZ() - exp.getMinMZ()));

    //----------------------------------------------------------------
    //Do the actual resampling
    BilinearInterpolation<double, double> bilip;
    bilip.getData().getEigenMatrix().resize(rows, cols);
    bilip.getData().getEigenMatrix().setZero();

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
    int scans = (int) bilip.getData().rows();
    int peaks = (int) bilip.getData().cols();

    bool use_log = getFlag_("log_intensity");

    MultiGradient gradient;
    String gradient_str = getStringOption_("gradient");
    if (!gradient_str.empty())
    {
      gradient.fromString(String("Linear|") + gradient_str);
    }
    else if (use_log)
    {
      gradient = MultiGradient::getDefaultGradientLogarithmicIntensityMode();
    }
    else
    {
      gradient = MultiGradient::getDefaultGradientLinearIntensityMode();
    }

    QImage image(peaks, scans, QImage::Format_RGB32);
    string s = getStringOption_("background_color");
    QColor background_color(s.c_str());

    string feature_color_string = getStringOption_("feature_color");
    QColor feature_color(feature_color_string.c_str());

    QPainter* painter = new QPainter(&image);
    painter->setPen(background_color);
    painter->fillRect(0, 0, peaks, scans, Qt::SolidPattern);
    delete painter;

    double factor = getDoubleOption_("max_intensity");
    if (factor == 0)
    {
      factor = bilip.getData().getEigenMatrix().maxCoeff();
    }
    // with a user-supplied gradient, we need to logarithmize explicitly;
    // by default, the gradient itself is adjusted to the log-scale:
    use_log &= !gradient_str.empty();
    if (use_log) factor = std::log(factor);

    factor /= 100.0;
    for (int i = 0; i < scans; ++i)
    {
      for (int j = 0; j < peaks; ++j)
      {
        double value = bilip.getData()(i, j);
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

    if (show_precursors)
    {
      markMS2Locations_(exp, image, getFlag_("transpose"),
                        getStringOption_("precursor_color").toQString(),
                        Size(getIntOption_("precursor_size")));
    }

    if (!in_featureXML.empty())
    {
      FeatureMap feature_map;
      FileHandler().loadFeatures(in_featureXML, feature_map, {FileTypes::FEATUREXML});
      markFeatureLocations_(feature_map, exp, image, getFlag_("transpose"), feature_color);
    }

    if (image.save(out.toQString(), format.c_str())) return EXECUTION_OK;
    else return CANNOT_WRITE_OUTPUT_FILE;
  }

};


int main(int argc, const char** argv)
{
  TOPPImageCreator tool;
  return tool.main(argc, argv);
}

/// @endcond
