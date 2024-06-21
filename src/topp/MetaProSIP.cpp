// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/PROCESSING/SCALING/Normalizer.h>
#include <OpenMS/FEATUREFINDER/GaussModel.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/SYSTEM/File.h>


#include <boost/math/distributions/normal.hpp>

#include <QtCore/QStringList>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtCore/QProcess>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#include <cmath>

//#define DEBUG_METAPROSIP

using namespace OpenMS;
using namespace std;
using boost::math::normal;

typedef map<double, double> MapRateToScoreType;
typedef pair<double, vector<double> > IsotopePattern;
typedef vector<IsotopePattern> IsotopePatterns;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MetaProSIP MetaProSIP

@brief Performs proteinSIP on peptide features for elemental flux analysis.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MetaProSIP.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MetaProSIP.html
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class MetaProSIP :
  public TOPPBase
{
public:
  MetaProSIP()
    : TOPPBase("MetaProSIP", "Performs proteinSIP on peptide features for elemental flux analysis."),
    ADDITIONAL_ISOTOPES(5),
    FEATURE_STRING("feature"),
    UNASSIGNED_ID_STRING("id"),
    UNIDENTIFIED_STRING("unidentified")
  {
  }

struct RateScorePair
{
  double rate = -1.;
  double score = -1.;
};

/// datastructure for reporting an incorporation event
struct SIPIncorporation
{
  double rate = -1.; ///< rate

  double correlation = -1.; ///< correlation coefficient

  double abundance = -1.; ///< abundance of isotopologue
#ifdef DEBUG_METAPROSIP
  PeakSpectrum theoretical; ///< peak spectrum as generated from the theoretical isotopic distribution. Large memory consumption.
#endif
};

/// datastructure for reporting a peptide with one or more incorporation rates
struct SIPPeptide
{
  AASequence sequence; ///< sequence of the peptide

  vector<String> accessions; ///< protein accessions of the peptide

  bool unique = true; ///< if the peptide is unique and therefor identifies the protein umambiguously

  double mz_theo = -1.; ///< theoretical mz

  double mass_theo = -1.; ///< uncharged theoretical mass

  double score = -1.; ///< search engine score or q-value if fdr filtering is applied

  double feature_rt = -1.; ///< measurement time of feature apex [s]

  double feature_mz = -1.; ///< mz of feature apex [s]

  //Size feature_scan_number; ///< scan number

  Int charge = 0; ///< charge of the peptide feature

  double mass_diff = 0.; // 13C or 15N mass difference

  double global_LR = -1.; ///< labeling ratio for the whole spectrum used to detect global drifts. 13C/(12C+13C) intensities. (15N analogous)

  vector<RateScorePair> correlation_maxima;

  MapRateToScoreType decomposition_map; // all rate to decomposition scores for the peptide

  MapRateToScoreType correlation_map; // all rate to correlation scores for the peptide

  double RR = -1.; ///< R squared of NNLS fit

  double explained_TIC_fraction = -1.; ///< fraction of the MS2 TIC that is explained by the maximum correlating decomposition weights

  String feature_type; ///< used to distinguish features from FeatureFinder, or synthetised from ids or averagine ids in reporting

  Size non_zero_decomposition_coefficients = 0; ///< decomposition coefficients significantly larger than 0

  PeakSpectrum reconstruction; ///< signal reconstruction (debugging)

  vector<double> reconstruction_monoistopic; ///< signal reconstruction of natural peptide (at mono-isotopic peak)

  PeakSpectrum accumulated;

  vector<SIPIncorporation> incorporations;

  IsotopePatterns patterns;

#ifdef DEBUG_METAPROSIP
  vector<PeakSpectrum> pattern_spectra;
#endif
};

///< comparator for vectors of SIPPeptides based on their size. Used to sort by group size.
struct SizeLess
{
  inline bool operator()(const vector<SIPPeptide>& a, const vector<SIPPeptide>& b) const
  {
    return a.size() < b.size();
  }

};

struct SequenceLess
{
  inline bool operator()(const pair<SIPPeptide, Size>& a, const pair<SIPPeptide, Size>& b) const
  {
    return a.first.sequence.toString() < b.first.sequence.toString();
  }

};

struct RIALess
{
  inline bool operator()(const SIPIncorporation& a, const SIPIncorporation& b) const
  {
    return a.rate < b.rate;
  }

};

class MetaProSIPInterpolation
{
public:
  ///< Determine score maxima from rate to score distribution using derivatives from spline interpolation
  static vector<RateScorePair> getHighPoints(double threshold, const MapRateToScoreType& rate2score, bool debug = false)
  {
    vector<RateScorePair> high_points;
    vector<double> x, y;

    // set proper boundaries (uniform spacing)
    x.push_back(-100.0 / (double)rate2score.size());
    y.push_back(0);

    // copy data
    for (MapRateToScoreType::const_iterator it = rate2score.begin(); it != rate2score.end(); ++it)
    {
      x.push_back(it->first);
      y.push_back(it->second);
    }

    if (rate2score.find(100.0) == rate2score.end() && x[x.size() - 1] < 100.0)
    {
      x.push_back(100.0);
      y.push_back(0);
    }

    const size_t n = x.size();

    //gte::IntpAkimaNonuniform1<double> spline(x.size(), &x.front(), &y.front());
    CubicSpline2d spline(x, y);

    if (debug)
    {
      OPENMS_LOG_DEBUG << x[0] << " " << x[n - 1] << " " << n << endl;
    }

    double last_dxdy = 0;
    for (double xi = x[0]; xi < x[n - 1]; xi += 0.01)
    {
      double dxdy = spline.derivatives(xi, 1);
      double yi = spline.eval(xi);

      if (debug)
      {
        cout << x[0] << " " << x[n - 1] << " " << xi << " " << yi << endl;
      }

      if (last_dxdy > 0.0 && dxdy <= 0 && yi > threshold)
      {
        RateScorePair rsp{};
        rsp.rate = xi;
        rsp.score = yi;
        high_points.push_back(rsp);
      }
      last_dxdy = dxdy;
    }

    if (debug)
    {
      OPENMS_LOG_DEBUG << "Found: " << high_points.size() << " local maxima." << endl;
      for (Size i = 0; i != high_points.size(); ++i)
      {
        OPENMS_LOG_DEBUG << high_points[i].rate << " " << high_points[i].score << endl;
      }
    }

    return high_points;
  }

};

class MetaProSIPClustering
{
public:
  static vector<double> getRIAClusterCenter(const vector<SIPPeptide>& sip_peptides, bool debug = false)
  {
    vector<double> cluster;
    MapRateToScoreType hist;

    for (vector<SIPPeptide>::const_iterator cit = sip_peptides.begin(); cit != sip_peptides.end(); ++cit)
    {
      // build histogram of rates
      for (vector<SIPIncorporation>::const_iterator iit = cit->incorporations.begin(); iit != cit->incorporations.end(); ++iit)
      {
        if (hist.find(iit->rate) == hist.end())
        {
          hist[iit->rate] = 1.0;
        }
        else
        {
          hist[iit->rate] += 1.0;
        }
      }
    }

    // kernel density estimation, TODO: binary search for 5 sigma boundaries
    vector<double> density(101, 0);
    for (Size i = 0; i != density.size(); ++i)
    {
      double sum = 0;
      for (MapRateToScoreType::const_iterator mit = hist.begin(); mit != hist.end(); ++mit)
      {
        normal s(mit->first, 2.0);
        sum += mit->second * pdf(s, (double)i);
      }
      density[i] = sum;
    }

    MapRateToScoreType ria_density;

    for (Size i = 0; i != density.size(); ++i)
    {
      ria_density[i] = density[i];
    }

    vector<RateScorePair> cluster_center = MetaProSIPInterpolation::getHighPoints(0.5, ria_density, debug);

    // return cluster centers
    for (vector<RateScorePair>::const_iterator cit = cluster_center.begin(); cit != cluster_center.end(); ++cit)
    {
      cluster.push_back(cit->rate);
    }
    return cluster;
  }

  //@Note sip peptides get reordered in same order as clusters
  static vector<vector<SIPPeptide> > clusterSIPPeptides(const vector<double>& centers, vector<SIPPeptide>& sip_peptides)
  {
    // one cluster for each cluster center
    vector<vector<SIPPeptide> > clusters(centers.size(), vector<SIPPeptide>());

    // assign sip peptide to cluster center with largest RIA
    for (vector<SIPPeptide>::const_iterator sit = sip_peptides.begin(); sit != sip_peptides.end(); ++sit)
    {
      const vector<SIPIncorporation>& incs = sit->incorporations;
      if (!incs.empty())
      {
        double largest_ria = incs[incs.size() - 1].rate;
        Size closest_cluster_idx = 0;
        double closest_cluster_dist = std::numeric_limits<double>::max();
        for (Size i = 0; i != centers.size(); ++i)
        {
          double dist = std::fabs(centers[i] - largest_ria);
          if (dist < closest_cluster_dist)
          {
            closest_cluster_dist = dist;
            closest_cluster_idx = i;
          }
        }

        // add SIP peptide to closest cluster
        clusters[closest_cluster_idx].push_back(*sit);
      }
    }

    // rearrange SIP peptides to reflect new order
    sip_peptides.clear();
    for (vector<vector<SIPPeptide> >::const_iterator sit = clusters.begin(); sit != clusters.end(); ++sit)
    {
      sip_peptides.insert(sip_peptides.end(), sit->begin(), sit->end());
    }

    return clusters;
  }

};

class MetaProSIPReporting
{
public:
  static void plotHeatMap(const String& output_dir, const String& tmp_path, const String& file_suffix, const String& file_extension, const vector<vector<double> >& binned_ria, vector<String> class_labels, Size debug_level = 0, const QString& executable = QString("R"))
  {
    String filename = String("heatmap") + file_suffix + "." + file_extension;
    String script_filename = String("heatmap") + file_suffix + String(".R");

    TextFile current_script;
    StringList ria_list, col_labels;

    for (Size i = 0; i != binned_ria[0].size(); ++i)
    {
      String label = String(i * (100 / binned_ria[0].size())) + "%-" + String((i + 1) * (100 / binned_ria[0].size())) + "%";
      col_labels.push_back(label);
    }

    for (vector<vector<double> >::const_iterator pit = binned_ria.begin(); pit != binned_ria.end(); ++pit)
    {
      for (vector<double>::const_iterator rit = pit->begin(); rit != pit->end(); ++rit)
      {
        ria_list.push_back(String(*rit));
      }
    }

    // row labels
    StringList row_labels;
    if (!class_labels.empty())
    {
      for (Size i = 0; i != class_labels.size(); ++i)
      {
        row_labels.push_back(class_labels[i]);
      }
    }

    // plot heatmap
    current_script.addLine("library(gplots)");
    String ria_list_string;
    ria_list_string.concatenate(ria_list.begin(), ria_list.end(), ",");
    current_script.addLine("mdat <- matrix(c(" + ria_list_string + "), ncol=" + String(binned_ria[0].size()) + ", byrow=TRUE)");

    if (file_extension == "png")
    {
      current_script.addLine("png('" + tmp_path + "/" + filename + "', width=1000, height=" + String(10 * binned_ria.size()) + ")");
    }
    else if (file_extension == "svg")
    {
      current_script.addLine("svg('" + tmp_path + "/" + filename + "', width=8, height=4.5)");
    }
    else if (file_extension == "pdf")
    {
      current_script.addLine("pdf('" + tmp_path + "/" + filename + "', width=8, height=4.5)");
    }

    String labRowString;
    if (row_labels.empty())
    {
      labRowString = "FALSE";
    }
    else
    {
      String row_labels_string;
      row_labels_string.concatenate(row_labels.begin(), row_labels.end(), "\",\"");
      labRowString = String("c(\"") + row_labels_string + "\")";
    }

    String col_labels_string;
    col_labels_string.concatenate(col_labels.begin(), col_labels.end(), "\",\"");

    current_script.addLine(R"(heatmap.2(mdat, dendrogram="none", col=colorRampPalette(c("black","red")), Rowv=FALSE, Colv=FALSE, key=FALSE, labRow=)" + labRowString + ",labCol=c(\"" + col_labels_string + R"("),trace="none", density.info="none"))");

    current_script.addLine("tmp<-dev.off()");
    current_script.store(tmp_path + "/" + script_filename);

    QProcess p;
    QStringList env = QProcess::systemEnvironment();
    env << QString("R_LIBS=") + tmp_path.toQString();
    p.setEnvironment(env);

    QStringList qparam;
    qparam << "--vanilla";
    if (debug_level < 1)
    {
      qparam << "--quiet";
    }
    qparam << "--slave" << "--file=" + QString(tmp_path.toQString() + "/" + script_filename.toQString());
    p.start(executable, qparam);
    p.waitForFinished(-1);
    int status = p.exitCode();

    // cleanup
    if (status != 0)
    {
      std::cerr << "Error: Process returned with non 0 status." << std::endl;
    }
    else
    {
      QFile(QString(tmp_path.toQString() + "/" + filename.toQString())).copy(output_dir.toQString() + "/heatmap" + file_suffix.toQString() + "." + file_extension.toQString());
      if (debug_level < 1)
      {
        QFile(QString(tmp_path.toQString() + "/" + script_filename.toQString())).remove();
        QFile(QString(tmp_path.toQString() + "/" + filename.toQString())).remove();
      }
    }
  }

  static void plotFilteredSpectra(const String& output_dir, const String& tmp_path, const String& file_suffix, const String& file_extension, const vector<SIPPeptide>& sip_peptides, Size debug_level = 0, const QString& executable = QString("R"))
  {
    String filename = String("spectrum_plot") + file_suffix + "." + file_extension;
    String script_filename = String("spectrum_plot") + file_suffix + String(".R");

    for (Size i = 0; i != sip_peptides.size(); ++i)
    {
      TextFile current_script;
      StringList mz_list;
      StringList intensity_list;

      for (Size j = 0; j != sip_peptides[i].accumulated.size(); ++j)
      {
        const Peak1D& peak = sip_peptides[i].accumulated[j];
        mz_list.push_back(String(peak.getMZ()));
        intensity_list.push_back(String(peak.getIntensity()));
      }

      String mz_list_string;
      mz_list_string.concatenate(mz_list.begin(), mz_list.end(), ",");

      String intensity_list_string;
      intensity_list_string.concatenate(intensity_list.begin(), intensity_list.end(), ",");

      current_script.addLine("mz<-c(" + mz_list_string + ")");
      current_script.addLine("int<-c(" + intensity_list_string + ")");
      current_script.addLine("x0=mz; x1=mz; y0=rep(0, length(x0)); y1=int");

      if (file_extension == "png")
      {
        current_script.addLine("png('" + tmp_path + "/" + filename + "')");
      }
      else if (file_extension == "svg")
      {
        current_script.addLine("svg('" + tmp_path + "/" + filename + "', width=8, height=4.5)");
      }
      else if (file_extension == "pdf")
      {
        current_script.addLine("pdf('" + tmp_path + "/" + filename + "', width=8, height=4.5)");
      }

      current_script.addLine("plot.new()");
      current_script.addLine("plot.window(xlim=c(min(mz),max(mz)), ylim=c(0,max(int)))");
      current_script.addLine("axis(1); axis(2)");
      current_script.addLine("title(xlab=\"m/z\")");
      current_script.addLine("title(ylab=\"intensity\")");
      current_script.addLine("box()");
      current_script.addLine("segments(x0,y0,x1,y1)");
      current_script.addLine("tmp<-dev.off()");
      current_script.store(tmp_path + "/" + script_filename);

      QProcess p;
      QStringList env = QProcess::systemEnvironment();
      env << QString("R_LIBS=") + tmp_path.toQString();
      p.setEnvironment(env);

      QStringList qparam;
      qparam << "--vanilla" << "--quiet" << "--slave" << "--file=" + QString(tmp_path.toQString() + "/" + script_filename.toQString());
      p.start(executable, qparam);
      p.waitForFinished(-1);
      int status = p.exitCode();

      if (status != 0)
      {
        std::cerr << "Error: Process returned with non 0 status." << std::endl;
      }
      else
      {
        QFile(QString(tmp_path.toQString() + "/" + filename.toQString())).copy(output_dir.toQString() + "/spectrum" + file_suffix.toQString() + "_rt_" + String(sip_peptides[i].feature_rt).toQString() + "." + file_extension.toQString());
        if (debug_level < 1)
        {
          QFile(QString(tmp_path.toQString() + "/" + script_filename.toQString())).remove();
          QFile(QString(tmp_path.toQString() + "/" + filename.toQString())).remove();
        }
      }
    }
  }

  static void writeHTML(const String& qc_output_directory, const String& file_suffix, const String& file_extension, const vector<SIPPeptide>& sip_peptides)
  {
    TextFile current_script;

    // html header
    current_script.addLine("<!DOCTYPE html>\n<html>\n<body>\n");

    // peptide heat map plot
    current_script.addLine(String("<h1>") + "peptide heat map</h1>");
    String peptide_heatmap_plot_filename = String("heatmap_peptide") + file_suffix + String(".") + file_extension;
    current_script.addLine("<p> <img src=\"" + peptide_heatmap_plot_filename + R"(" alt="graphic"></p>)");

    for (Size i = 0; i != sip_peptides.size(); ++i)
    {
      // heading
      current_script.addLine(String("<h1>") + "RT: " + String(sip_peptides[i].feature_rt) + "</h1>");

      current_script.addLine("<table border=\"1\">");
      // sequence table row
      current_script.addLine("<tr>");
      current_script.addLine("<td>sequence</td>");
      current_script.addLine(String("<td>") + sip_peptides[i].sequence.toString() + "</td>");
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>rt (min.)</td>");
      current_script.addLine(String("<td>" + String::number(sip_peptides[i].feature_rt / 60.0, 2) + "</td>"));
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>rt (sec.)</td>");
      current_script.addLine(String("<td>" + String::number(sip_peptides[i].feature_rt, 2) + "</td>"));
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>mz</td>");
      current_script.addLine(String("<td>" + String::number(sip_peptides[i].feature_mz, 4) + "</td>"));
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>theo. mz</td>");
      current_script.addLine(String("<td>" + String::number(sip_peptides[i].mz_theo, 4) + "</td>"));
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>charge</td>");
      current_script.addLine(String("<td>" + String(sip_peptides[i].charge) + "</td>"));
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>feature type</td>");
      current_script.addLine(String("<td>" + String(sip_peptides[i].feature_type) + "</td>"));
      current_script.addLine("</tr>");

      if (!sip_peptides[i].accessions.empty())
      {
        current_script.addLine(String("<tr>"));
        current_script.addLine("<td>accessions</td>");
        current_script.addLine(String("<td>" + *sip_peptides[i].accessions.begin() + "</td>"));
        current_script.addLine(String("</tr>"));

        current_script.addLine(String("<tr>"));
        current_script.addLine("<td>unique</td>");
        current_script.addLine(String("<td>" + String(sip_peptides[i].unique) + "</td>"));
        current_script.addLine(String("</tr>"));
      }

      current_script.addLine(String("<tr>"));
      current_script.addLine("<td>search score</td>");
      current_script.addLine(String("<td>") + String(sip_peptides[i].score) + "</td>");
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>global labeling ratio</td>");
      current_script.addLine(String("<td>") + String::number(sip_peptides[i].global_LR, 2) + "</td>");
      current_script.addLine("</tr>");

      current_script.addLine("<tr>");
      current_script.addLine("<td>R squared</td>");
      current_script.addLine(String("<td>") + String::number(sip_peptides[i].RR, 2) + "</td>");
      current_script.addLine("</tr>");

      current_script.addLine("</table>");

      // table header of incorporations
      current_script.addLine("<p>");
      current_script.addLine("<table border=\"1\">");
      current_script.addLine("<tr>");
      for (Size k = 0; k != sip_peptides[i].incorporations.size(); ++k)
      {
        current_script.addLine(String("<td>RIA") + String(k + 1) + "</td>");
        current_script.addLine(String("<td>CORR.") + String(k + 1) + "</td>");
        current_script.addLine(String("<td>INT") + String(k + 1) + "</td>");
      }
      current_script.addLine("</tr>");

      // table of incorporations
      current_script.addLine("<tr>");
      for (Size k = 0; k != sip_peptides[i].incorporations.size(); ++k)
      {
        SIPIncorporation p = sip_peptides[i].incorporations[k];
        current_script.addLine(String("<td>") + String::number(p.rate, 2) + "</td>");
        current_script.addLine(String("<td>") + String::number(p.correlation, 2) + "</td>");
        current_script.addLine(String("<td>") + String::number(p.abundance, 0) + "</td>");
      }
      current_script.addLine("</tr>");

      current_script.addLine("</table>");

      // spectrum plot
      String spectrum_filename = String("spectrum") + file_suffix + "_rt_" + String(sip_peptides[i].feature_rt) + "." + file_extension;
      current_script.addLine("<p> <img src=\"" + spectrum_filename + R"(" alt="graphic"></p>)");

      // score plot
      String score_filename = String("scores") + file_suffix + "_rt_" + String(sip_peptides[i].feature_rt) + "." + file_extension;
      current_script.addLine("<p> <img src=\"" + score_filename + R"(" alt="graphic"></p>)");
    }
    current_script.addLine("\n</body>\n</html>");
    current_script.store(qc_output_directory.toQString() + "/index" + file_suffix.toQString() + ".html");
  }

  static void plotScoresAndWeights(const String& output_dir, const String& tmp_path, const String& file_suffix, const String& file_extension, const vector<SIPPeptide>& sip_peptides, double score_plot_yaxis_min, Size debug_level = 0, const QString& executable = QString("R"))
  {
    String score_filename = String("score_plot") + file_suffix + file_extension;
    String script_filename = String("score_plot") + file_suffix + String(".R");

    for (Size i = 0; i != sip_peptides.size(); ++i)
    {
      TextFile current_script;
      StringList rate_dec_list;
      StringList rate_corr_list;
      StringList weights_list;
      StringList corr_list;

      for (MapRateToScoreType::const_iterator mit = sip_peptides[i].decomposition_map.begin(); mit != sip_peptides[i].decomposition_map.end(); ++mit)
      {
        rate_dec_list.push_back(String(mit->first));
        weights_list.push_back(String(mit->second));
      }

      for (MapRateToScoreType::const_iterator mit = sip_peptides[i].correlation_map.begin(); mit != sip_peptides[i].correlation_map.end(); ++mit)
      {
        rate_corr_list.push_back(String(mit->first));
        corr_list.push_back(String(mit->second));
      }

      String rate_dec_list_string;
      rate_dec_list_string.concatenate(rate_dec_list.begin(), rate_dec_list.end(), ",");

      String weights_list_string;
      weights_list_string.concatenate(weights_list.begin(), weights_list.end(), ",");

      String rate_corr_list_string;
      rate_corr_list_string.concatenate(rate_corr_list.begin(), rate_corr_list.end(), ",");

      String corr_list_string;
      corr_list_string.concatenate(corr_list.begin(), corr_list.end(), ",");

      current_script.addLine("rate_dec<-c(" + rate_dec_list_string + ")");
      current_script.addLine("dec<-c(" + weights_list_string + ")");
      current_script.addLine("if (max(dec)!=0) {dec<-dec/max(dec)}");
      current_script.addLine("rate_corr<-c(" + rate_corr_list_string + ")");
      current_script.addLine("corr<-c(" + corr_list_string + ")");

      if (score_plot_yaxis_min >= 0)
      {
        current_script.addLine("corr[corr<0]=0"); // truncate at 0 for better drawing
      }

      current_script.addLine("x0=rate_dec; x1=rate_dec; y0=rep(0, length(x0)); y1=dec"); // create R segments for decomposition score (vertical bars)
      if (file_extension == "png")
      {
        current_script.addLine("png('" + tmp_path + "/" + score_filename + "')");
      }
      else if (file_extension == "svg")
      {
        current_script.addLine("svg('" + tmp_path + "/" + score_filename + "', width=8, height=4.5)");
      }
      else if (file_extension == "pdf")
      {
        current_script.addLine("pdf('" + tmp_path + "/" + score_filename + "', width=8, height=4.5)");
      }
      current_script.addLine("plot.new()");
      current_script.addLine("plot.window(xlim=c(0,100), ylim=c(" + String(score_plot_yaxis_min) + ",1))");
      current_script.addLine("axis(1); axis(2)");
      current_script.addLine("title(xlab=\"RIA\")");
      current_script.addLine("title(ylab=\"normalized weight / corr.\")");
      current_script.addLine("box()");
      current_script.addLine("segments(x0,y0,x1,y1, col='red')");
      current_script.addLine("lines(x=rate_corr, y=corr, col='blue')");
      current_script.addLine("legend('bottomright', horiz=FALSE, xpd=TRUE, col=c('red', 'blue'), lwd=2, c('weights', 'correlation'))");
      current_script.addLine("tmp<-dev.off()");
      current_script.store(tmp_path + "/" + script_filename);

      QProcess p;
      QStringList env = QProcess::systemEnvironment();
      env << QString("R_LIBS=") + tmp_path.toQString();
      p.setEnvironment(env);

      QStringList qparam;
      qparam << "--vanilla" << "--quiet" << "--slave" << "--file=" + QString(tmp_path.toQString() + "/" + script_filename.toQString());
      p.start(executable, qparam);
      p.waitForFinished(-1);
      int status = p.exitCode();

      if (status != 0)
      {
        std::cerr << "Error: Process returned with non 0 status." << std::endl;
      }
      else
      {
        QFile(QString(tmp_path.toQString() + "/" + score_filename.toQString())).copy(output_dir.toQString() + "/scores" + file_suffix.toQString() + "_rt_" + String(sip_peptides[i].feature_rt).toQString() + "." + file_extension.toQString());
        if (debug_level < 1)
        {
          QFile(QString(tmp_path.toQString() + "/" + script_filename.toQString())).remove();
          QFile(QString(tmp_path.toQString() + "/" + score_filename.toQString())).remove();
        }
      }
    }
  }

  static void createQualityReport(const String& tmp_path,
                                  const String& qc_output_directory,
                                  const String& file_suffix,
                                  const String& file_extension,
                                  const vector<vector<SIPPeptide> >& sip_peptide_cluster,
                                  Size n_heatmap_bins,
                                  double score_plot_y_axis_min,
                                  bool report_natural_peptides,
                                  const QString& executable = QString("R"))
  {
    vector<SIPPeptide> sip_peptides;
    for (vector<vector<SIPPeptide> >::const_iterator cit = sip_peptide_cluster.begin(); cit != sip_peptide_cluster.end(); ++cit)
    {
      for (vector<SIPPeptide>::const_iterator sit = cit->begin(); sit != cit->end(); ++sit)
      {
        // skip non natural peptides for reporting if flag is set
        if (!report_natural_peptides && sit->incorporations.size() == 1 && sit->incorporations[0].rate < 5.0)
        {
          continue;
        }
        sip_peptides.push_back(*sit);
      }
    }

    // heat map based on peptide RIAs
    OPENMS_LOG_INFO << "Plotting peptide heat map of " << sip_peptides.size() << endl;
    vector<vector<double> > binned_peptide_ria;
    vector<String> class_labels;
    createBinnedPeptideRIAData_(n_heatmap_bins, sip_peptide_cluster, binned_peptide_ria, class_labels);
    plotHeatMap(qc_output_directory, tmp_path, "_peptide" + file_suffix, file_extension, binned_peptide_ria, class_labels, 0, executable);

    OPENMS_LOG_INFO << "Plotting filtered spectra for quality report" << endl;
    plotFilteredSpectra(qc_output_directory, tmp_path, file_suffix, file_extension, sip_peptides, 0, executable);

    OPENMS_LOG_INFO << "Plotting correlation score and weight distribution" << endl;
    plotScoresAndWeights(qc_output_directory, tmp_path, file_suffix, file_extension, sip_peptides, score_plot_y_axis_min, 0, executable);

    if (file_extension != "pdf") // html doesn't support pdf as image
    {
      writeHTML(qc_output_directory, file_suffix, file_extension, sip_peptides);
    }
  }

  static void createCSVReport(vector<vector<SIPPeptide> >& sippeptide_cluster, ofstream& os, map<String, String>& proteinid_to_description)
  {
    SVOutStream out_csv_stream(os, "\t", "_", String::NONE);
    // sort clusters by non increasing size
    sort(sippeptide_cluster.rbegin(), sippeptide_cluster.rend(), SizeLess());

    for (Size pep_clust_i = 0; pep_clust_i != sippeptide_cluster.size(); ++pep_clust_i)
    {
      const vector<SIPPeptide>& current_cluster = sippeptide_cluster[pep_clust_i];

      // Group
      map<String, vector<SIPPeptide> > all_peptides; // map sequence to SIPPeptide
      map<String, vector<SIPPeptide> > ambigous_peptides; // map sequence to SIPPeptide
      map<String, map<String, vector<SIPPeptide> > > unambigous_proteins; // map Accession to unmodified String to SIPPeptides

      for (Size k = 0; k != current_cluster.size(); ++k)
      {
        const SIPPeptide& current_SIPpeptide = current_cluster[k];
        String seq = current_SIPpeptide.sequence.toUnmodifiedString();
        if (current_SIPpeptide.unique)
        {
          String first_accession = *current_SIPpeptide.accessions.begin();
          unambigous_proteins[first_accession][seq].push_back(current_SIPpeptide);
        }
        else
        {
          ambigous_peptides[current_SIPpeptide.sequence.toUnmodifiedString()].push_back(current_SIPpeptide);
        }
        all_peptides[seq].push_back(current_SIPpeptide);
      }

      Size n_all_peptides = all_peptides.size(); // # of different (on sequence level) unique and non-unique peptides
      //Size n_ambigous_peptides = ambigous_peptides.size();
      Size n_unambigous_proteins = unambigous_proteins.size();

      // determine median global LR of whole group
      vector<double> group_global_LRs;
      vector<double> group_number_RIAs;
      for (map<String, vector<SIPPeptide> >::const_iterator all_it = all_peptides.begin(); all_it != all_peptides.end(); ++all_it)
      {
        for (vector<SIPPeptide>::const_iterator v_it = all_it->second.begin(); v_it != all_it->second.end(); ++v_it)
        {
          group_global_LRs.push_back(v_it->global_LR);
          group_number_RIAs.push_back(v_it->incorporations.size());
        }
      }
      double group_global_LR = Math::median(group_global_LRs.begin(), group_global_LRs.end(), false);

      Size group_number_RIA = lround(Math::median(group_number_RIAs.begin(), group_number_RIAs.end(), false)); // median number of RIAs
      // Group header
      // Distinct peptides := different (on sequence level) unique and non-unique peptides
      out_csv_stream << String("Group ") + String(pep_clust_i + 1) << "# Distinct Peptides" << "# Unambiguous Proteins" << "Median Global LR";
      for (Size i = 0; i != group_number_RIA; ++i)
      {
        out_csv_stream << "median RIA " + String(i + 1);
      }
      out_csv_stream << endl;

      out_csv_stream << "" << n_all_peptides << n_unambigous_proteins << group_global_LR;

      // collect 1th, 2nd, ... RIA of the group based on the peptide RIAs
      vector<vector<double> > group_RIAs(group_number_RIA, vector<double>());
      vector<double> group_RIA_medians(group_number_RIA, 0);

      for (map<String, vector<SIPPeptide> >::const_iterator all_it = all_peptides.begin(); all_it != all_peptides.end(); ++all_it)
      {
        for (vector<SIPPeptide>::const_iterator v_it = all_it->second.begin(); v_it != all_it->second.end(); ++v_it)
        {
          for (Size i = 0; i != group_number_RIA; ++i)
          {
            if (i == v_it->incorporations.size())
            {
              break;
            }
            group_RIAs[i].push_back(v_it->incorporations[i].rate);
          }
        }
      }

      for (Size i = 0; i != group_number_RIA; ++i)
      {
        group_RIA_medians[i] = Math::median(group_RIAs[i].begin(), group_RIAs[i].end(), false);
      }

      for (Size i = 0; i != group_number_RIA; ++i)
      {
        out_csv_stream << String(group_RIA_medians[i]);
      }
      out_csv_stream << endl;

      // unambiguous protein level
      for (map<String, map<String, vector<SIPPeptide> > >::const_iterator prot_it = unambigous_proteins.begin(); prot_it != unambigous_proteins.end(); ++prot_it)
      {
        // determine median global LR of protein
        vector<double> protein_global_LRs;
        vector<double> protein_number_RIAs;
        for (map<String, vector<SIPPeptide> >::const_iterator pept_it = prot_it->second.begin(); pept_it != prot_it->second.end(); ++pept_it)
        {
          for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
          {
            protein_global_LRs.push_back(v_it->global_LR);
            protein_number_RIAs.push_back(v_it->incorporations.size());
          }
        }
        double protein_global_LR = Math::median(protein_global_LRs.begin(), protein_global_LRs.end(), false);
        Size protein_number_RIA = (Size)(Math::median(protein_number_RIAs.begin(), protein_number_RIAs.end(), false) + 0.5); // median number of RIAs

        out_csv_stream << "" << "Protein Accession" << "Description" << "# Unique Peptides" << "Median Global LR";
        for (Size i = 0; i != protein_number_RIA; ++i)
        {
          out_csv_stream << "median RIA " + String(i + 1);
        }
        out_csv_stream << endl;

        String protein_accession = prot_it->first;
        String protein_description = "none";
        if (proteinid_to_description.find(protein_accession.trim().toUpper()) != proteinid_to_description.end())
        {
          protein_description = proteinid_to_description.at(protein_accession.trim().toUpper());
        }

        out_csv_stream << "" << protein_accession << protein_description << prot_it->second.size() << protein_global_LR;

        vector<vector<double> > protein_RIAs(protein_number_RIA, vector<double>());
        vector<double> protein_RIA_medians(protein_number_RIA, 0);

        // ratio to natural decomposition
        vector<vector<double> > protein_ratio(protein_number_RIA, vector<double>());
        vector<double> protein_ratio_medians(protein_number_RIA, 0);

        // collect 1th, 2nd, ... RIA of the protein based on the peptide RIAs
        for (map<String, vector<SIPPeptide> >::const_iterator pept_it = prot_it->second.begin(); pept_it != prot_it->second.end(); ++pept_it)
        {
          for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
          {
            for (Size i = 0; i != protein_number_RIA; ++i)
            {
              if (i == v_it->incorporations.size())
              {
                break;
              }
              protein_RIAs[i].push_back(v_it->incorporations[i].rate);
              protein_ratio[i].push_back(v_it->incorporations[i].abundance);
            }
          }
        }

        for (Size i = 0; i != protein_number_RIA; ++i)
        {
          protein_RIA_medians[i] = Math::median(protein_RIAs[i].begin(), protein_RIAs[i].end(), false);
          protein_ratio_medians[i] = Math::median(protein_ratio[i].begin(), protein_ratio[i].end(), false);
        }

        for (Size i = 0; i != protein_number_RIA; ++i)
        {
          out_csv_stream << String(protein_RIA_medians[i]);
        }

        out_csv_stream << endl;

        // print header of unique peptides
        out_csv_stream << "" << "" << "Peptide Sequence" << "RT" << "Exp. m/z" << "Theo. m/z" << "Charge" << "Score" << "TIC fraction" << "#non-natural weights" << "";
        Size max_incorporations = 0;
        for (map<String, vector<SIPPeptide> >::const_iterator pept_it = prot_it->second.begin(); pept_it != prot_it->second.end(); ++pept_it)
        {
          for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
          {
            max_incorporations = std::max(v_it->incorporations.size(), max_incorporations);
          }
        }

        for (Size i = 0; i != max_incorporations; ++i)
        {
          out_csv_stream << "RIA " + String(i + 1) << "INT " + String(i + 1) << "Cor. " + String(i + 1);
        }
        out_csv_stream << "Peak intensities" << "Global LR" << endl;

        // print data of unique peptides
        for (map<String, vector<SIPPeptide> >::const_iterator pept_it = prot_it->second.begin(); pept_it != prot_it->second.end(); ++pept_it)
        {
          for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
          {
            out_csv_stream << "" << "" << v_it->sequence.toString() << String::number(v_it->feature_rt / 60.0, 2) << String::number(v_it->feature_mz, 4) << v_it->mz_theo << v_it->charge << v_it->score << v_it->explained_TIC_fraction << v_it->non_zero_decomposition_coefficients << "";
            for (vector<SIPIncorporation>::const_iterator incorps = v_it->incorporations.begin(); incorps != v_it->incorporations.end(); ++incorps)
            {
              out_csv_stream << String::number(incorps->rate, 1) << String::number(incorps->abundance, 0) << String::number(incorps->correlation, 2);
            }

            // blank entries for nicer formatting
            for (Int q = 0; q < (Int)max_incorporations - (Int)v_it->incorporations.size(); ++q)
            {
              out_csv_stream << "" << "" << "";
            }

            // output peak intensities
            String peak_intensities;
            for (PeakSpectrum::const_iterator p = v_it->accumulated.begin(); p != v_it->accumulated.end(); ++p)
            {
              peak_intensities += String::number(p->getIntensity(), 0) + " ";
            }
            out_csv_stream << peak_intensities;
            out_csv_stream << v_it->global_LR;

            out_csv_stream << endl;
          }
        }
      }

      // print header of non-unique peptides below the protein section
      Size max_incorporations = 0;
      for (map<String, vector<SIPPeptide> >::const_iterator pept_it = ambigous_peptides.begin(); pept_it != ambigous_peptides.end(); ++pept_it)
      {
        for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
        {
          max_incorporations = std::max(v_it->incorporations.size(), max_incorporations);
        }
      }

      out_csv_stream << "Non-Unique Peptides" << "Accessions" << "Peptide Sequence" << "Descriptions" << "Score" << "RT" << "Exp. m/z" << "Theo. m/z" << "Charge" << "#non-natural weights" << "";

      for (Size m = 0; m != max_incorporations; ++m)
      {
        out_csv_stream << "RIA " + String(m + 1) << "INT " + String(m + 1) << "Cor. " + String(m + 1);
      }
      out_csv_stream << "Peak intensities" << "Global LR" << endl;

      // print data of non-unique peptides below the protein section
      for (map<String, vector<SIPPeptide> >::const_iterator pept_it = ambigous_peptides.begin(); pept_it != ambigous_peptides.end(); ++pept_it)
      {
        // build up the protein accession string for non-unique peptides. Only the first 3 accessions are added.
        for (vector<SIPPeptide>::const_iterator v_it = pept_it->second.begin(); v_it != pept_it->second.end(); ++v_it)
        {
          String accessions_string;
          String description_string = "none";

          for (Size ac = 0; ac != v_it->accessions.size(); ++ac)
          {
            if (ac >= 3) // only print at most 3 accessions as these can be quite numorous
            {
              accessions_string += "...";
              break;
            }
            String protein_accession = v_it->accessions[ac];
            accessions_string += protein_accession;

            if (proteinid_to_description.find(protein_accession.trim().toUpper()) != proteinid_to_description.end())
            {
              if (description_string == "none")
              {
                description_string = "";
              }
              description_string += proteinid_to_description.at(protein_accession.trim().toUpper());
            }

            if (ac < v_it->accessions.size() - 1)
            {
              accessions_string += ", ";
              if (description_string != "none")
              {
                description_string += ", ";
              }
            }
          }

          out_csv_stream << "" << accessions_string << v_it->sequence.toString() << description_string << v_it->score << String::number(v_it->feature_rt / 60.0, 2) << String::number(v_it->feature_mz, 4) << v_it->mz_theo << v_it->charge << v_it->non_zero_decomposition_coefficients << "";

          // output variable sized RIA part
          for (vector<SIPIncorporation>::const_iterator incorps = v_it->incorporations.begin(); incorps != v_it->incorporations.end(); ++incorps)
          {
            out_csv_stream << String::number(incorps->rate, 1) << String::number(incorps->abundance, 0) << String::number(incorps->correlation, 2);
          }

          // blank entries for nicer formatting
          for (Int q = 0; q < (Int)max_incorporations - (Int)v_it->incorporations.size(); ++q)
          {
            out_csv_stream << "" << "" << "";
          }

          // output peak intensities
          String peak_intensities;
          for (PeakSpectrum::const_iterator p = v_it->accumulated.begin(); p != v_it->accumulated.end(); ++p)
          {
            peak_intensities += String::number(p->getIntensity(), 0) + " ";
          }
          out_csv_stream << peak_intensities;
          out_csv_stream << v_it->global_LR;
          out_csv_stream << endl;
        }
      }
    }
    os.close();
  }

  static void createPeptideCentricCSVReport(const String& in_mzML, const String& file_extension, vector<vector<SIPPeptide> >& sippeptide_cluster, ofstream& os, map<String, String>& proteinid_to_description, String qc_output_directory, String file_suffix, bool report_natural_peptides)
  {
    SVOutStream out_csv_stream(os, "\t", "_", String::NONE);

    // sort clusters by non increasing size
    sort(sippeptide_cluster.rbegin(), sippeptide_cluster.rend(), SizeLess());

    // store SIP peptide with cluster index for peptide centric view on data
    vector<pair<SIPPeptide, Size> > peptide_to_cluster_index;
    for (Size i = 0; i != sippeptide_cluster.size(); ++i)
    {
      const vector<SIPPeptide>& current_cluster = sippeptide_cluster[i];
      for (Size k = 0; k != current_cluster.size(); ++k)
      {
        peptide_to_cluster_index.emplace_back(current_cluster[k], i);
      }
    }

    OPENMS_LOG_INFO << "Writing " << peptide_to_cluster_index.size() << " peptides to peptide centric csv." << endl;

    // sort by sequence
    sort(peptide_to_cluster_index.begin(), peptide_to_cluster_index.end(), SequenceLess());

    out_csv_stream << "Peptide Sequence" << "Feature" << "Quality Report Spectrum" << "Quality report scores" << "Sample Name" << "Protein Accessions" << "Description" << "Unique" << "#Ambiguity members"
                   << "Score" << "RT" << "Exp. m/z" << "Theo. m/z" << "Charge" << "TIC fraction" << "#non-natural weights" << "Peak intensities" << "Group" << "Global Peptide LR";

    for (Size i = 1; i <= 10; ++i)
    {
      out_csv_stream << "RIA " + String(i) << "LR of RIA " + String(i) << "INT " + String(i) << "Cor. " + String(i);
    }
    out_csv_stream << std::endl;

    for (Size i = 0; i != peptide_to_cluster_index.size(); ++i)
    {
      const SIPPeptide& current_SIPpeptide = peptide_to_cluster_index[i].first;

      // skip non natural peptides for repoting if flag is set
      if (!report_natural_peptides
        && current_SIPpeptide.incorporations.size() == 1
        && current_SIPpeptide.incorporations[0].rate < 5.0)
      {
        continue;
      }

      const Size& current_cluster_index = peptide_to_cluster_index[i].second;

      // output peptide sequence
      out_csv_stream << current_SIPpeptide.sequence.toString() << current_SIPpeptide.feature_type;

      // output quality report links if available
      if (qc_output_directory.empty() || file_suffix.empty()) // if no qc plots have been generated or no unique file_suffix has been provided we can't generate links to spectra and scores
      {
        out_csv_stream << "" << "" << in_mzML;
      }
      else
      {
        String qr_spectrum_filename = String("file://") + qc_output_directory + "/" + String("spectrum") + file_suffix + "_rt_" + String(current_SIPpeptide.feature_rt) + "." + file_extension;
        String qr_scores_filename = String("file://") + qc_output_directory + "/" + String("scores") + file_suffix + "_rt_" + String(current_SIPpeptide.feature_rt) + "." + file_extension;
        out_csv_stream << qr_spectrum_filename << qr_scores_filename << in_mzML;
      }

      // output protein accessions and descriptions
      String accession_string;
      String protein_descriptions = "none";
      for (Size j = 0; j != current_SIPpeptide.accessions.size(); ++j)
      {
        String current_accession = current_SIPpeptide.accessions[j];
        current_accession.trim().toUpper();
        accession_string += current_accession;

        if (proteinid_to_description.find(current_accession) != proteinid_to_description.end())
        {
          if (protein_descriptions == "none")
          {
            protein_descriptions = proteinid_to_description.at(current_accession);
          }
          else
          {
            protein_descriptions += proteinid_to_description.at(current_accession);
          }
        }

        // add "," between accessions
        if (j != current_SIPpeptide.accessions.size() - 1)
        {
          accession_string += ",";
          protein_descriptions += ",";
        }
      }

      out_csv_stream << accession_string << protein_descriptions << current_SIPpeptide.unique << current_SIPpeptide.accessions.size() << current_SIPpeptide.score << String::number(current_SIPpeptide.feature_rt / 60.0, 2)
                     << String::number(current_SIPpeptide.feature_mz, 4) << String::number(current_SIPpeptide.mz_theo, 4) << current_SIPpeptide.charge << current_SIPpeptide.explained_TIC_fraction << current_SIPpeptide.non_zero_decomposition_coefficients;

      // output peak intensities
      String peak_intensities;
      for (PeakSpectrum::const_iterator p = current_SIPpeptide.accumulated.begin(); p != current_SIPpeptide.accumulated.end(); ++p)
      {
        peak_intensities += String::number(p->getIntensity(), 0) + " ";
      }
      out_csv_stream << peak_intensities;

      out_csv_stream << current_cluster_index << current_SIPpeptide.global_LR;

      for (Size j = 0; j != current_SIPpeptide.incorporations.size(); ++j)
      {
        const double ria = current_SIPpeptide.incorporations[j].rate;
        const double abundance = current_SIPpeptide.incorporations[j].abundance;
        const double corr = current_SIPpeptide.incorporations[j].correlation;

        double LR_of_RIA = 0;
        if (ria < 1.5) // first RIA hast natural abundance
        {
          LR_of_RIA = abundance / current_SIPpeptide.incorporations[0].abundance;
        }
        out_csv_stream << String::number(ria, 1) << String::number(LR_of_RIA, 1) << String::number(abundance, 1) << String::number(corr, 1);
      }
      out_csv_stream << endl;
    }

    out_csv_stream << endl;
    os.close();
  }

protected:
  static void createBinnedPeptideRIAData_(const Size n_heatmap_bins, const vector<vector<SIPPeptide> >& sip_clusters, vector<vector<double> >& binned_peptide_ria, vector<String>& cluster_labels)
  {
    cluster_labels.clear();
    binned_peptide_ria.clear();

    for (vector<vector<SIPPeptide> >::const_iterator cit = sip_clusters.begin(); cit != sip_clusters.end(); ++cit)
    {
      const vector<SIPPeptide>& sip_peptides = *cit;
      for (vector<SIPPeptide>::const_iterator pit = sip_peptides.begin(); pit != sip_peptides.end(); ++pit)
      {
        vector<double> binned(n_heatmap_bins, 0.0);
        for (vector<SIPIncorporation>::const_iterator iit = pit->incorporations.begin(); iit != pit->incorporations.end(); ++iit)
        {
          Int bin = static_cast<Int>(iit->rate / 100.0 * n_heatmap_bins);
          bin = bin > (Int)binned.size() - 1 ? (Int)binned.size() - 1 : bin;
          bin = bin < 0 ? 0 : bin;
          binned[bin] = log1p(iit->abundance);
        }
        binned_peptide_ria.push_back(binned);
        cluster_labels.push_back((String)(cit - sip_clusters.begin()));
      }
    }
  }

};

class MetaProSIPDecomposition
{
public:
  ///> Perform the decomposition
  static Int calculateDecompositionWeightsIsotopicPatterns(Size n_bins, const vector<double>& isotopic_intensities, const IsotopePatterns& patterns, MapRateToScoreType& map_rate_to_decomposition_weight, SIPPeptide& sip_peptide)
  {
    Matrix<double> beta(n_bins, 1);
    Matrix<double> intensity_vector(isotopic_intensities.size(), 1);

    for (Size p = 0; p != isotopic_intensities.size(); ++p)
    {
      intensity_vector(p, 0) = isotopic_intensities[p];
    }

    Matrix<double> basis_matrix(isotopic_intensities.size(), n_bins);

    for (Size row = 0; row != isotopic_intensities.size(); ++row)
    {
      for (Size col = 0; col != n_bins; ++col)
      {
        const vector<double>& pattern = patterns[col].second;
        if (row <= n_bins)
        {
          basis_matrix(row, col) = pattern[row];
        }
        else
        {
          basis_matrix(row, col) = 0;
        }
      }
    }

    Int result = NonNegativeLeastSquaresSolver::solve(basis_matrix, intensity_vector, beta);

    for (Size p = 0; p != n_bins; ++p)
    {
      map_rate_to_decomposition_weight[(double)p / n_bins * 100.0] = beta(p, 0);
    }

    // calculate R squared
    double S_tot = 0;
    double mean = accumulate(isotopic_intensities.begin(), isotopic_intensities.end(), 0.0) / isotopic_intensities.size();
    for (Size row = 0; row != isotopic_intensities.size(); ++row)
    {
      S_tot += pow(isotopic_intensities[row] - mean, 2);
    }

    double S_err = 0;
    PeakSpectrum reconstructed;

    for (Size row = 0; row != isotopic_intensities.size(); ++row)
    {
      double predicted = 0;
      for (Size col = 0; col != n_bins; ++col)
      {
        predicted += basis_matrix(row, col) * beta(col, 0);
      }
      Peak1D peak;
      peak.setIntensity(predicted);
      peak.setMZ(sip_peptide.mz_theo + sip_peptide.mass_diff / sip_peptide.charge * row);
      reconstructed.push_back(peak);
      S_err += pow(isotopic_intensities[row] - predicted, 2);
    }

    for (Size row = 0; row != 5; ++row)
    {
      double predicted = 0;
      for (Size col = 0; col != 3; ++col)
      {
        predicted += basis_matrix(row, col) * beta(col, 0);
      }
      sip_peptide.reconstruction_monoistopic.push_back(predicted);
    }

    sip_peptide.RR = 1.0 - (S_err / S_tot);
    sip_peptide.reconstruction = reconstructed;

    return result;
  }

  // Template calculations for base matrix

  ///> Given a peptide sequence calculate the theoretical isotopic patterns given all incorporations rate (13C Version)
  ///> extend isotopic patterns by additional_isotopes to collect other element higher isotopes at 100% incorporation
  static IsotopePatterns calculateIsotopePatternsFor13CRange(const AASequence& peptide, Size additional_isotopes = 5)
  {
    IsotopePatterns ret;
    const Element* e1 = ElementDB::getInstance()->getElement("Carbon");
    Element* e2 = const_cast<Element*>(e1);

    EmpiricalFormula peptide_ef = peptide.getFormula();
    Size MAXISOTOPES = static_cast<Size>(peptide_ef.getNumberOf(e1));

    // calculate empirical formula of modifications - these can not be labeled via substrate feeding and must be taken care of in pattern calculation
    AASequence unmodified_peptide = AASequence::fromString(peptide.toUnmodifiedString());
    EmpiricalFormula unmodified_peptide_ef = unmodified_peptide.getFormula();
    UInt max_labeling_carbon = (UInt)unmodified_peptide_ef.getNumberOf(e1); // max. number of atoms that can be labeled
    EmpiricalFormula modifications_ef = peptide_ef - unmodified_peptide_ef; // difference formula for modifications (note that it can contain positive/negative numbers)

    if (modifications_ef.getNumberOf(e1) > 0) // modification adds additional (unlabeled) carbon atoms
    {
      IsotopeDistribution modification_dist = modifications_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_carbon + additional_isotopes));

      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)max_labeling_carbon)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(12, 1.0 - a);
        isotopes.insert(13, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = unmodified_peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_carbon + additional_isotopes));
        dist.set(CoarseIsotopePatternGenerator().convolve(dist.getContainer(), modification_dist.getContainer())); // convolve with modification distribution (which follows the natural distribution)
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }
    else
    {

      // calculate isotope distribution for a given peptide and varying incorporation rates
      // modification of isotope distribution in static ElementDB
      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)MAXISOTOPES)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(12, 1.0 - a);
        isotopes.insert(13, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(MAXISOTOPES + additional_isotopes));

        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;

        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(12, 0.9893f);
    isotopes.insert(13, 0.0107f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static Size getNumberOfLabelingElements(const String& labeling_element, const AASequence& peptide)
  {
    const Element * e;
    if (labeling_element == "N")
    {
      e = ElementDB::getInstance()->getElement("Nitrogen");
    }
    else if (labeling_element == "C")
    {
      e = ElementDB::getInstance()->getElement("Carbon");
    }
    else if (labeling_element == "H")
    {
      e = ElementDB::getInstance()->getElement("Hydrogen");
    }
    else if (labeling_element == "O")
    {
      e = ElementDB::getInstance()->getElement("Oxygen");
    }
    else
    {
      return 0;
    }

    // try to determine if modification adds or removes elements
    AASequence unmodified_peptide = AASequence::fromString(peptide.toUnmodifiedString());
    EmpiricalFormula unmodified_peptide_ef = unmodified_peptide.getFormula();
    int labeling_element_mods_excluded = unmodified_peptide_ef.getNumberOf(e);

    EmpiricalFormula peptide_ef = peptide.getFormula();
    int labeling_element_mods_included = peptide_ef.getNumberOf(e);

    int diff = labeling_element_mods_included - labeling_element_mods_excluded;

    if (diff >= 0) // common case, mod added unlabeled elements
    {
      return labeling_element_mods_excluded;
    }
    else // special case, mod results in loss of labeling element
    {
      return labeling_element_mods_included;
    }
  }

  ///> Given a peptide sequence calculate the theoretical isotopic patterns given all incorporations rate (15C Version)
  ///> extend isotopic patterns by additional_isotopes to collect other element higher isotopes at 100% incorporation
  static IsotopePatterns calculateIsotopePatternsFor15NRange(const AASequence& peptide, Size additional_isotopes = 5)
  {
    IsotopePatterns ret;

    const Element* e1 = ElementDB::getInstance()->getElement("Nitrogen");
    Element* e2 = const_cast<Element*>(e1);

    EmpiricalFormula peptide_ef = peptide.getFormula();
    UInt MAXISOTOPES = static_cast<UInt>(peptide_ef.getNumberOf(e1));

    // calculate empirical formula of modifications - these can not be labeled via substrate feeding and must be taken care of in pattern calculation
    AASequence unmodified_peptide = AASequence::fromString(peptide.toUnmodifiedString());
    EmpiricalFormula unmodified_peptide_ef = unmodified_peptide.getFormula();
    UInt max_labeling_nitrogens = (UInt)unmodified_peptide_ef.getNumberOf(e1); // max. number of nitrogen atoms that can be labeled
    EmpiricalFormula modifications_ef = peptide_ef - unmodified_peptide_ef; // difference formula for modifications (note that it can contain positive/negative numbers)

    if (modifications_ef.getNumberOf(e1) > 0) // modification adds additional (unlabeled) nitrogen atoms
    {
      IsotopeDistribution modification_dist = modifications_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_nitrogens + additional_isotopes));
      for (double abundance = 0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)max_labeling_nitrogens)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(14, 1.0 - a);
        isotopes.insert(15, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = unmodified_peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_nitrogens + additional_isotopes));
        dist.set(CoarseIsotopePatternGenerator().convolve(dist.getContainer(), modification_dist.getContainer())); // calculate convolution with isotope distribution of modification(s)
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }
    else
    {
      // calculate isotope distribution for a given peptide and varying incoperation rates
      // modification of isotope distribution in static ElementDB
      for (double abundance = 0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)MAXISOTOPES)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(14, 1.0 - a);
        isotopes.insert(15, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(MAXISOTOPES + additional_isotopes));
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }
    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(14, 0.99632f);
    isotopes.insert(15, 0.368f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor2HRange(const AASequence& peptide, Size additional_isotopes = 5)
  {
    IsotopePatterns ret;

    const Element* e1 = ElementDB::getInstance()->getElement("Hydrogen");
    Element* e2 = const_cast<Element*>(e1);

    EmpiricalFormula peptide_ef = peptide.getFormula();
    Size MAXISOTOPES = static_cast<Size>(peptide_ef.getNumberOf(e1));

    // calculate empirical formula of modifications - these can not be labeled via substrate feeding and must be taken care of in pattern calculation
    AASequence unmodified_peptide = AASequence::fromString(peptide.toUnmodifiedString());
    EmpiricalFormula unmodified_peptide_ef = unmodified_peptide.getFormula();
    UInt max_labeling_element = (UInt)unmodified_peptide_ef.getNumberOf(e1); // max. number of atoms that can be labeled
    EmpiricalFormula modifications_ef = peptide_ef - unmodified_peptide_ef; // difference formula for modifications (note that it can contain positive/negative numbers)

    if (modifications_ef.getNumberOf(e1) > 0) // modification adds additional (unlabeled) atoms
    {
      IsotopeDistribution modification_dist = modifications_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_element + additional_isotopes));
      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)max_labeling_element)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(1, 1.0 - a);
        isotopes.insert(2, a);
        e2->setIsotopeDistribution(isotopes);

        IsotopeDistribution dist = unmodified_peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_element + additional_isotopes));
        dist.set(CoarseIsotopePatternGenerator().convolve(dist.getContainer(), modification_dist.getContainer())); // convole with modification distribution (which follows the natural distribution)
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }
    else
    {
      // calculate isotope distribution for a given peptide and varying incoperation rates
      // modification of isotope distribution in static ElementDB
      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)MAXISOTOPES)
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(1, 1.0 - a);
        isotopes.insert(2, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(MAXISOTOPES + additional_isotopes));
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(1, 0.999885f);
    isotopes.insert(2, 0.000115f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor18ORange(const AASequence& peptide, Size additional_isotopes = 5)
  {
    IsotopePatterns ret;

    const Element* e1 = ElementDB::getInstance()->getElement("Oxygen");
    Element* e2 = const_cast<Element*>(e1);

    EmpiricalFormula peptide_ef = peptide.getFormula();
    Size MAXISOTOPES = static_cast<Size>(peptide_ef.getNumberOf(e1));
    // calculate empirical formula of modifications - these can not be labeled via substrate feeding and must be taken care of in pattern calculation
    AASequence unmodified_peptide = AASequence::fromString(peptide.toUnmodifiedString());
    EmpiricalFormula unmodified_peptide_ef = unmodified_peptide.getFormula();
    UInt max_labeling_element = (UInt)unmodified_peptide_ef.getNumberOf(e1); // max. number of atoms that can be labeled
    EmpiricalFormula modifications_ef = peptide_ef - unmodified_peptide_ef; // difference formula for modifications (note that it can contain positive/negative numbers)

    if (modifications_ef.getNumberOf(e1) > 0) // modification adds additional (unlabeled) atoms
    {
      IsotopeDistribution modification_dist = modifications_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_element + additional_isotopes));
      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / static_cast<double>(max_labeling_element * 2.0))
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.insert(1, 1.0 - a);
        isotopes.insert(2, 0.0); // 17O is neglectable (=0.038%)
        isotopes.insert(3, a);
        e2->setIsotopeDistribution(isotopes);

        IsotopeDistribution dist = unmodified_peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_labeling_element * 2 + additional_isotopes)); // 2 * isotopic traces
        dist.set(CoarseIsotopePatternGenerator().convolve(dist.getContainer(), modification_dist.getContainer())); // convole with modification distribution (which follows the natural distribution)
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }
    else
    {
      // calculate isotope distribution for a given peptide and varying incoperation rates
      // modification of isotope distribution in static ElementDB
      for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / static_cast<double>(MAXISOTOPES * 2.0))
      {
        double a = abundance / 100.0;
        IsotopeDistribution isotopes;
        isotopes.clear();
        isotopes.insert(1, 1.0 - a);
        isotopes.insert(2, 0.0); // 17O is neglectable (=0.038%)
        isotopes.insert(3, a);
        e2->setIsotopeDistribution(isotopes);
        IsotopeDistribution dist = peptide_ef.getIsotopeDistribution(CoarseIsotopePatternGenerator(MAXISOTOPES * 2 + additional_isotopes)); // 2 * isotopic traces
        IsotopeDistribution::ContainerType container = dist.getContainer();
        vector<double> intensities;
        for (Size i = 0; i != container.size(); ++i)
        {
          intensities.push_back(container[i].getIntensity());
        }
        ret.push_back(make_pair(abundance, intensities));
      }
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(1, 0.99757f);
    isotopes.insert(2, 0.00038f);
    isotopes.insert(3, 0.00205f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor15NRangeOfAveraginePeptide(double mass)
  {
    IsotopePatterns ret;
    const Element* e1 = ElementDB::getInstance()->getElement("Nitrogen");
    Element* e2 = const_cast<Element*>(e1);

    // calculate number of expected labeling elements using averagine model
    Size element_count = static_cast<Size>(mass * 0.0122177302837372);

    // calculate isotope distribution for a given peptide and varying incoperation rates
    // modification of isotope distribution in static ElementDB
    for (double abundance = 0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)element_count)
    {
      double a = abundance / 100.0;
      IsotopeDistribution isotopes;
      isotopes.clear();
      isotopes.insert(14, 1.0 - a);
      isotopes.insert(15, a);
      e2->setIsotopeDistribution(isotopes);
      CoarseIsotopePatternGenerator solver(element_count);
      auto dist = solver.estimateFromPeptideWeight(mass);
      IsotopeDistribution::ContainerType container = dist.getContainer();
      vector<double> intensities;
      for (Size i = 0; i != container.size(); ++i)
      {
        intensities.push_back(container[i].getIntensity());
      }
      ret.push_back(make_pair(abundance, intensities));
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(14, 0.99632f);
    isotopes.insert(15, 0.368f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor13CRangeOfAveraginePeptide(double mass)
  {
    IsotopePatterns ret;
    const Element* e1 = ElementDB::getInstance()->getElement("Carbon");
    Element* e2 = const_cast<Element*>(e1);
    Size element_count = static_cast<Size>(mass * 0.0444398894906044);

    // calculate isotope distribution for a given peptide and varying incoperation rates
    // modification of isotope distribution in static ElementDB
    for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)element_count)
    {
      double a = abundance / 100.0;
      IsotopeDistribution isotopes;
      isotopes.clear();
      isotopes.insert(12, 1.0 - a);
      isotopes.insert(13, a);
      e2->setIsotopeDistribution(isotopes);
      CoarseIsotopePatternGenerator solver(element_count);
      auto dist = solver.estimateFromPeptideWeight(mass);
      IsotopeDistribution::ContainerType container = dist.getContainer();
      vector<double> intensities;
      for (Size i = 0; i != container.size(); ++i)
      {
        intensities.push_back(container[i].getIntensity());
      }
      ret.push_back(make_pair(abundance, intensities));
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.insert(12, 0.9893f);
    isotopes.insert(13, 0.010f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor2HRangeOfAveraginePeptide(double mass)
  {
    IsotopePatterns ret;

    const Element* e1 = ElementDB::getInstance()->getElement("Hydrogen");
    Element* e2 = const_cast<Element*>(e1);
    Size element_count = static_cast<Size>(mass * 0.06981572169);

    // calculate isotope distribution for a given peptide and varying incoperation rates
    // modification of isotope distribution in static ElementDB
    for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)element_count)
    {
      double a = abundance / 100.0;
      IsotopeDistribution isotopes;
      isotopes.clear();
      isotopes.insert(1, 1.0 - a);
      isotopes.insert(2, a);
      e2->setIsotopeDistribution(isotopes);
      CoarseIsotopePatternGenerator solver(element_count);
      auto dist = solver.estimateFromPeptideWeight(mass);
      IsotopeDistribution::ContainerType container = dist.getContainer();
      vector<double> intensities;
      for (Size i = 0; i != container.size(); ++i)
      {
        intensities.push_back(container[i].getIntensity());
      }
      ret.push_back(make_pair(abundance, intensities));
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(1, 0.999885f);
    isotopes.insert(2, 0.000115f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }

  static IsotopePatterns calculateIsotopePatternsFor18ORangeOfAveraginePeptide(double mass)
  {
    IsotopePatterns ret;

    const Element* e1 = ElementDB::getInstance()->getElement("Oxygen");
    Element* e2 = const_cast<Element*>(e1);
    Size element_count = static_cast<Size>(mass * 0.01329399039);

    // calculate isotope distribution for a given peptide and varying incoperation rates
    // modification of isotope distribution in static ElementDB
    for (double abundance = 0.0; abundance < 100.0 - 1e-8; abundance += 100.0 / (double)element_count)
    {
      double a = abundance / 100.0;
      IsotopeDistribution isotopes;
      isotopes.clear();
      isotopes.insert(1, 1.0 - a);
      isotopes.insert(2, 0);
      isotopes.insert(3, a);
      e2->setIsotopeDistribution(isotopes);
      CoarseIsotopePatternGenerator solver(element_count * 2); // spaces are 2 Da between 18O and 16O but we observe isotopic peaks at every (approx.) nominal mass
      auto dist = solver.estimateFromPeptideWeight(mass);
      IsotopeDistribution::ContainerType container = dist.getContainer();
      vector<double> intensities;
      for (Size i = 0; i != container.size(); ++i)
      {
        intensities.push_back(container[i].getIntensity());
      }
      ret.push_back(make_pair(abundance, intensities));
    }

    // reset to natural occurance
    IsotopeDistribution isotopes;
    isotopes.clear();
    isotopes.insert(1, 0.99757f);
    isotopes.insert(2, 0.00038f);
    isotopes.insert(3, 0.00205f);
    e2->setIsotopeDistribution(isotopes);
    return ret;
  }
};

class MetaProSIPXICExtraction
{
public:
  static vector<vector<double> > extractXICs(double seed_rt, vector<double> xic_mzs, double mz_toelrance_ppm, double rt_tolerance_s, const PeakMap& peak_map)
  {
    // point on first spectrum in tolerance window
    PeakMap::ConstIterator rt_begin = peak_map.RTBegin(seed_rt - rt_tolerance_s);

    // point on after last spectrum in tolerance window
    PeakMap::ConstIterator rt_end = peak_map.RTBegin(seed_rt + rt_tolerance_s);

    // create set containing all rts of spectra in tolerance window
    set<double> all_rts;
    for (PeakMap::ConstIterator rt_it = rt_begin; rt_it != rt_end; ++rt_it)
    {
      all_rts.insert(rt_it->getRT());
    }

    vector<vector<double> > xics(xic_mzs.size(), vector<double>());

    for (Size i = 0; i < xic_mzs.size(); ++i)
    {
      // create and initialize xic to contain values for all rts
      map<double, double> xic; // rt to summed intensity
      for (set<double>::const_iterator sit = all_rts.begin(); sit != all_rts.end(); ++sit)
      {
        xic[*sit] = 0;
      }

      double mz_da = mz_toelrance_ppm * xic_mzs[i] * 1e-6; // mz tolerance in Dalton
      PeakMap::ConstAreaIterator it = peak_map.areaBeginConst(seed_rt - rt_tolerance_s, seed_rt + rt_tolerance_s, xic_mzs[i] - mz_da, xic_mzs[i] + mz_da);

      for (; it != peak_map.areaEndConst(); ++it)
      {
        double rt = it.getRT();
        if (xic.find(rt) != xic.end())
        {
          xic[rt] += it->getIntensity();
        }
        else
        {
          OPENMS_LOG_WARN << "RT: " << rt << " not contained in rt set." << endl;
        }
      }

      // copy map to vector for easier processing
      vector<double> v;
      for (map<double, double>::const_iterator xic_it = xic.begin(); xic_it != xic.end(); ++xic_it)
      {
        v.push_back(xic_it->second);
      }

      xics[i] = v;
    }
    return xics;
  }

  static vector<double> correlateXICsToMono(const vector<vector<double> >& xics)
  {
    vector<double> rrs(xics.size(), 0); // correlation of isotopic xics to monoisotopic xic

    rrs[0] = 1.0; // perfect correlation of monoisotopic trace to itself

    for (Size i = 1; i < xics.size(); ++i)
    {
      rrs[i] = Math::pearsonCorrelationCoefficient(xics[0].begin(), xics[0].end(), xics[i].begin(), xics[i].end());
    }
    return rrs;
  }

  static vector<double> extractXICsOfIsotopeTraces(Size element_count, double mass_diff, double mz_tolerance_ppm, double rt_tolerance_s, double seed_rt, double seed_mz, double seed_charge, const PeakMap& peak_map, const double min_corr_mono = -1.0)
  {
    vector<double> xic_mzs;

    // calculate centers of XICs to be extracted
    for (Size k = 0; k != element_count; ++k)
    {
      double mz = seed_mz + k * mass_diff / seed_charge;
      xic_mzs.push_back(mz);
    }

    // extract xics
    vector<vector<double> > xics = extractXICs(seed_rt, xic_mzs, mz_tolerance_ppm, rt_tolerance_s, peak_map);

    vector<double> xic_intensities(xics.size(), 0.0);
    if (min_corr_mono > 0)
    {
      // calculate correlation to mono-isotopic peak
      vector<double> RRs = correlateXICsToMono(xics);

      // sum over XICs to yield one intensity value for each XIC. If correlation to mono-isotopic is lower then threshold, delete intensity.
      for (Size i = 0; i != xic_intensities.size(); ++i)
      {
        double v = std::accumulate(xics[i].begin(), xics[i].end(), 0.0);
        xic_intensities[i] = RRs[i] > min_corr_mono ? v : 0.0;
      }
    }
    else // correlation disabled so just take the XIC intensities
    {
      for (Size i = 0; i != xic_intensities.size(); ++i)
      {
        xic_intensities[i] = std::accumulate(xics[i].begin(), xics[i].end(), 0.0);
      }
    }

    return xic_intensities;
  }

};

class RIntegration
{
public:
  // Perform a simple check if R and all R dependencies are there
  static bool checkRDependencies(const String& tmp_path, StringList package_names, const QString& executable = QString("R"))
  {
    String random_name = String::random(8);
    String script_filename = tmp_path + String("/") + random_name + String(".R");

    // check if R in path and can be executed
    TextFile checkRInPath;
    checkRInPath.addLine("q()");
    checkRInPath.store(script_filename);

    OPENMS_LOG_INFO << "Checking R...";
    {
      QProcess p;
      p.setProcessChannelMode(QProcess::MergedChannels);
      QStringList env = QProcess::systemEnvironment();
      env << QString("R_LIBS=") + tmp_path.toQString();
      p.setEnvironment(env);

      QStringList checkRinPathQParam;
      checkRinPathQParam << "--vanilla" << "--quiet" << "--slave" << "--file=" + script_filename.toQString();
      p.start(executable, checkRinPathQParam);
      p.waitForFinished(-1);

      if (p.error() == QProcess::FailedToStart || p.exitStatus() == QProcess::CrashExit || p.exitCode() != 0)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Can't execute R. Do you have R installed? Check if the path to R is in your system path variable." << std::endl;
        return false;
      }
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    // check dependencies
    OPENMS_LOG_INFO << "Checking R dependencies. If package is not found we will try to install it in your temp directory...";
    TextFile current_script;
    current_script.addLine("LoadOrInstallPackage <-function(x)");
    current_script.addLine("{");
    current_script.addLine("  x <-as.character(substitute(x))");
    current_script.addLine("  if (isTRUE(x %in%.packages(all.available = TRUE)))");
    current_script.addLine("  {");
    current_script.addLine("    eval(parse(text = paste(\"library(\", x, \")\", sep = \"\")))");
    current_script.addLine("  }");
    current_script.addLine("  else");
    current_script.addLine("  {");
    current_script.addLine("    options(repos = structure(c(CRAN = \"http://cran.rstudio.com/\")))");
    current_script.addLine("    update.packages()");
    current_script.addLine("    eval(parse(text = paste(\"install.packages('\", x, \"')\", sep = \"\")))");
    current_script.addLine("    eval(parse(text = paste(\"library(\", x, \")\", sep = \"\")))");
    current_script.addLine("  }");
    current_script.addLine("}");
    for (StringList::const_iterator it = package_names.begin(); it != package_names.end(); ++it)
    {
      current_script.addLine("LoadOrInstallPackage(" + *it + ")");
    }

    current_script.store(script_filename);

    QProcess p;
    p.setProcessChannelMode(QProcess::MergedChannels);
    QStringList env = QProcess::systemEnvironment();
    env << QString("R_LIBS=") + tmp_path.toQString();
    p.setEnvironment(env);

    QStringList qparam;
    qparam << "--vanilla" << "--quiet" << "--slave" << "--file=" + script_filename.toQString();
    p.start(executable, qparam);
    p.waitForFinished(-1);
    int status = p.exitCode();

    if (status != 0)
    {
      OPENMS_LOG_ERROR << "\nProblem finding all R dependencies. Check if R and following libraries are installed:" << std::endl;
      for (TextFile::ConstIterator line_it = current_script.begin(); line_it != current_script.end(); ++line_it)
      {
        OPENMS_LOG_ERROR << *line_it  << std::endl;
      }
      QString s = p.readAllStandardOutput();
      OPENMS_LOG_ERROR << s.toStdString() << std::endl;
      return false;
    }
    OPENMS_LOG_INFO << " success" << std::endl;
    return true;
  }

};


protected:
  Size ADDITIONAL_ISOTOPES;
  std::string FEATURE_STRING;
  std::string UNASSIGNED_ID_STRING;
  std::string UNIDENTIFIED_STRING;
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in_mzML", "<file>", "", "Centroided MS1 data");
    setValidFormats_("in_mzML", ListUtils::create<String>("mzML"));

    registerInputFile_("in_fasta", "<file>", "", "Protein sequence database");
    setValidFormats_("in_fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_csv", "<file>", "", "Column separated file with feature fitting result.");
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));

    registerOutputFile_("out_peptide_centric_csv", "<file>", "", "Column separated file with peptide centric result.");
    setValidFormats_("out_peptide_centric_csv", ListUtils::create<String>("csv"));

    registerInputFile_("in_featureXML", "<file>", "", "Feature data annotated with identifications (IDMapper)");
    setValidFormats_("in_featureXML", ListUtils::create<String>("featureXML"));

    registerInputFile_("r_executable", "<file>", "R", "Path to the R executable (default: 'R')", false, false, {"is_executable"});

    registerDoubleOption_("mz_tolerance_ppm", "<tol>", 10.0, "Tolerance in ppm", false, true);

    registerDoubleOption_("rt_tolerance_s", "<tol>", 30.0, "Tolerance window around feature rt for XIC extraction", false, true);

    registerDoubleOption_("intensity_threshold", "<tol>", 10.0, "Intensity threshold to collect peaks in the MS1 spectrum.", false, true);

    registerDoubleOption_("correlation_threshold", "<tol>", 0.7, "Correlation threshold for reporting a RIA", false, true);

    registerDoubleOption_("xic_threshold", "<tol>", 0.7, "Minimum correlation to mono-isotopic peak for retaining a higher isotopic peak. If featureXML from reference file is used it should be disabled (set to -1) as no mono-isotopic peak is expected to be present.", false, true);

    registerDoubleOption_("decomposition_threshold", "<tol>", 0.7, "Minimum R-squared of decomposition that must be achieved for a peptide to be reported.", false, true);

    registerDoubleOption_("weight_merge_window", "<tol>", 5.0, "Decomposition coefficients within +- this rate window will be combined", false, true);

    registerDoubleOption_("min_correlation_distance_to_averagine", "<tol>", -1.0, "Minimum difference in correlation between incorporation pattern and averagine pattern. Positive values filter all RIAs passing the correlation threshold but that also show a better correlation to an averagine peptide. Disabled for values <= -1", false, true);

    registerDoubleOption_("pattern_15N_TIC_threshold", "<threshold>", 0.95, "The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.", false, true);
    registerDoubleOption_("pattern_13C_TIC_threshold", "<threshold>", 0.95, "The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.", false, true);
    registerDoubleOption_("pattern_2H_TIC_threshold", "<threshold>", 0.95, "The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.", false, true);
    registerDoubleOption_("pattern_18O_TIC_threshold", "<threshold>", 0.95, "The most intense peaks of the theoretical pattern contributing to at least this TIC fraction are taken into account.", false, true);
    registerIntOption_("heatmap_bins", "<threshold>", 20, "Number of RIA bins for heat map generation.", false, true);

    registerStringOption_("plot_extension", "<extension>", "png", "Extension used for plots (png|svg|pdf).", false, true);
    StringList valid_extensions;
    valid_extensions.push_back("png");
    valid_extensions.push_back("svg");
    valid_extensions.push_back("pdf");
    setValidStrings_("plot_extension", valid_extensions);

    registerStringOption_("qc_output_directory", "<directory>", "", "Output directory for the quality report", false, true);

    registerStringOption_("labeling_element", "<parameter>", "C", "Which element (single letter code) is labeled.", false);
    StringList valid_element;
    valid_element.push_back("C");
    valid_element.push_back("N");
    valid_element.push_back("H");
    valid_element.push_back("O");
    setValidStrings_("labeling_element", valid_element);

    registerFlag_("use_unassigned_ids", "Include identifications not assigned to a feature in pattern detection.", true);

    registerFlag_("use_averagine_ids", "Use averagine peptides as model to perform pattern detection on unidentified peptides.", true);

    registerFlag_("report_natural_peptides", "Whether purely natural peptides are reported in the quality report.", true);

    registerFlag_("filter_monoisotopic", "Try to filter out mono-isotopic patterns to improve detection of low RIA patterns", true);

    registerFlag_("cluster", "Perform grouping", true);

    registerDoubleOption_("observed_peak_fraction", "<threshold>", 0.5, "Fraction of observed/expected peaks.", false, true);

    registerIntOption_("min_consecutive_isotopes", "<threshold>", 2, "Minimum number of consecutive isotopic intensities needed.", false, true);

    registerDoubleOption_("score_plot_yaxis_min", "<threshold>", 0.0, "The minimum value of the score axis. Values smaller than zero usually only make sense if the observed peak fraction is set to 0.", false, true);

    registerStringOption_("collect_method", "<method>", "correlation_maximum", "How RIAs are collected.", false, true);
    StringList valid_collect_method;
    valid_collect_method.push_back("correlation_maximum");
    valid_collect_method.push_back("decomposition_maximum");
    setValidStrings_("collect_method", valid_collect_method);

    registerDoubleOption_("lowRIA_correlation_threshold", "<tol>", -1, "Correlation threshold for reporting low RIA patterns. Disable and take correlation_threshold value for negative values.", false, true);
  }

  ///> filter intensity to remove noise or additional incorporation peaks that otherwise might interfere with correlation calculation
  void filterIsotopicIntensities(vector<double>::const_iterator& pattern_begin, vector<double>::const_iterator& pattern_end,
                                 vector<double>::const_iterator& intensities_begin, vector<double>::const_iterator& intensities_end, double TIC_threshold = 0.99)
  {
    if (std::distance(pattern_begin, pattern_end) != std::distance(intensities_begin, intensities_end))
    {
      OPENMS_LOG_ERROR << "Error: size of pattern and collected intensities don't match!: (pattern " << std::distance(pattern_begin, pattern_end) << ") (intensities " << std::distance(intensities_begin, intensities_end) << ")" << endl;
    }

    if (pattern_begin == pattern_end)
    {
      return;
    }

    // determine order of peaks based on intensities
    vector<double>::const_iterator b_it = pattern_begin;
    vector<double>::const_iterator e_it = pattern_end;
    // create intensity to offset map for sorting
    vector<std::pair<double, Int> > intensity_to_offset;
    for (; b_it != e_it; ++b_it)
    {
      std::pair<double, Int> intensity_offset_pair = make_pair(*b_it, std::distance(pattern_begin, b_it));
      intensity_to_offset.push_back(intensity_offset_pair); // pair: intensity, offset to pattern_begin iterator
    }
    // sort by intensity (highest first)
    std::sort(intensity_to_offset.begin(), intensity_to_offset.end(), std::greater<pair<double, Int> >());

    // determine sequence of (neighbouring) peaks needed to achieve threshold * 100 % TIC in the patterns
    double TIC = 0.0;
    Int min_offset = std::distance(pattern_begin, pattern_end);
    Int max_offset = 0;

    for (vector<std::pair<double, Int> >::const_iterator it = intensity_to_offset.begin(); it != intensity_to_offset.end(); ++it)
    {
      TIC += it->first;
      if (it->second < min_offset)
      {
        min_offset = it->second;
      }

      if (it->second > max_offset)
      {
        max_offset = it->second;
      }

      if (TIC > TIC_threshold)
      {
        break;
      }
    }

    vector<double>::const_iterator tmp_pattern_it(pattern_begin);
    vector<double>::const_iterator tmp_intensity_it(intensities_begin);

    std::advance(pattern_begin, min_offset);
    std::advance(intensities_begin, min_offset);
    std::advance(tmp_pattern_it, max_offset + 1);
    std::advance(tmp_intensity_it, max_offset + 1);

    pattern_end = tmp_pattern_it;
    intensities_end = tmp_intensity_it;

    //cout << "after: " << std::distance(pattern_begin, pattern_end) << " " << min_offset << " " << max_offset << endl;
  }

  ///< Calculates the correlation between measured isotopic_intensities and the theoretical isotopic patterns for all incorporation rates
  void calculateCorrelation(Size n_element, const vector<double>& isotopic_intensities, IsotopePatterns patterns,
                            MapRateToScoreType& map_rate_to_correlation_score, String labeling_element, double mass, double min_correlation_distance_to_averagine)
  {
    double min_observed_peak_fraction = getDoubleOption_("observed_peak_fraction");

    if (debug_level_ > 0)
    {
      cout << "Calculating " << patterns.size() << " isotope patterns with " << ADDITIONAL_ISOTOPES << " additional isotopes." << endl;
    }

    double TIC_threshold(0.0);

    // N15 has smaller RIA resolution and multiple RIA peaks tend to overlap more in correlation. This reduces the width of the pattern leading to better distinction
    if (labeling_element == "N")
    {
      TIC_threshold = getDoubleOption_("pattern_15N_TIC_threshold");
    }
    else if (labeling_element == "C")
    {
      TIC_threshold = getDoubleOption_("pattern_13C_TIC_threshold");
    }
    else if (labeling_element == "H")
    {
      TIC_threshold = getDoubleOption_("pattern_2H_TIC_threshold");
    }
    else if (labeling_element == "O")
    {
      TIC_threshold = getDoubleOption_("pattern_18O_TIC_threshold");
    }
    double max_incorporation_rate = 100.0;
    double incorporation_step = max_incorporation_rate / (double)n_element;

    // calculate correlation with a natural averagine peptide (used to filter out coeluting peptides)
    double peptide_weight = mass;

    const Size AVERAGINE_CORR_OFFSET = 3;

    // calculate correlation for averagine peptides
    std::vector<double> averagine_correlation(isotopic_intensities.size(), 0.0);

    // extended by zeros on both sides to simplify correlation
    vector<double> ext_isotopic_intensities(AVERAGINE_CORR_OFFSET, 0.0);
    ext_isotopic_intensities.insert(ext_isotopic_intensities.end(), isotopic_intensities.begin(), isotopic_intensities.end());
    for (Size i = 0; i != AVERAGINE_CORR_OFFSET; ++i)
    {
      ext_isotopic_intensities.push_back(0.0);
    }

    for (Size ii = 0; ii < isotopic_intensities.size(); ++ii)
    {
      // calculate isotope distribution of averagine peptide as this will be used to detect spurious correlations with coeluting peptides
      // Note: actually it would be more accurate to use 15N-14N or 13C-12C distances. This doesn't affect averagine distribution much so this approximation is sufficient. (see TODO)
      double current_weight = peptide_weight + ii * 1.0; // TODO: use 13C-12C or 15N-14N instead of 1.0 as mass distance to be super accurate
      CoarseIsotopePatternGenerator solver(10);
      IsotopeDistribution averagine = solver.estimateFromPeptideWeight(current_weight);

      IsotopeDistribution::ContainerType averagine_intensities_pairs = averagine.getContainer();

      // zeros to the left for sliding window correlation
      std::vector<double> averagine_intensities(AVERAGINE_CORR_OFFSET, 0.0); // add 0 intensity bins left to actual averagine pattern

      for (Size i = 0; i != averagine_intensities_pairs.size(); ++i)
      {
        averagine_intensities.push_back(averagine_intensities_pairs[i].getIntensity());
      }

      // zeros to the right
      for (Size i = 0; i != AVERAGINE_CORR_OFFSET; ++i)
      {
        averagine_intensities.push_back(0.0);
      }

      // number of bins that can be correlated
      Int max_correlated_values = std::min((int)ext_isotopic_intensities.size() - ii, averagine_intensities.size());

      double corr_with_averagine = Math::pearsonCorrelationCoefficient(averagine_intensities.begin(), averagine_intensities.begin() + max_correlated_values,
                                                                       ext_isotopic_intensities.begin() + ii, ext_isotopic_intensities.begin() + ii + max_correlated_values);
      averagine_correlation[ii] = corr_with_averagine;
    }

    // calculate correlation of RIA peptide with measured data
    for (Size ii = 0; ii != patterns.size(); ++ii)
    {
      double rate = (double)ii * incorporation_step;

      vector<double>::const_iterator pattern_begin = patterns[ii].second.begin();
      vector<double>::const_iterator pattern_end = patterns[ii].second.end();
      vector<double>::const_iterator intensities_begin = isotopic_intensities.begin();
      vector<double>::const_iterator intensities_end = isotopic_intensities.end();

      filterIsotopicIntensities(pattern_begin, pattern_end, intensities_begin, intensities_end, TIC_threshold);
      Size zeros = 0;
      for (vector<double>::const_iterator it = intensities_begin; it != intensities_end; ++it)
      {
        if (*it < 1e-8)
        {
          zeros++;
        }
      }

      // remove correlations with only very few peaks
      if ((double)zeros / (double)std::distance(intensities_begin, intensities_end) > min_observed_peak_fraction)
      {
        map_rate_to_correlation_score[rate] = 0;
        continue;
      }

      double correlation_score = Math::pearsonCorrelationCoefficient(pattern_begin, pattern_end, intensities_begin, intensities_end);

      // remove correlations that show higher similarity to an averagine peptide
      if (rate > 5.0 && correlation_score < averagine_correlation[ii] + min_correlation_distance_to_averagine)
      {
        map_rate_to_correlation_score[rate] = 0;
        continue;
      }

      // cout << ii << "\t" << std::distance(intensities_end, intensities_begin) << "\t" << std::distance(intensities_begin, isotopic_intensities.begin()) << "\t" << std::distance(intensities_end, isotopic_intensities.begin()) << endl;
      if (std::isnan(correlation_score))
      {
        correlation_score = 0.0;
      }
      map_rate_to_correlation_score[rate] = correlation_score;
    }
  }

  ///< Returns highest scoring rate and score pair in the map
  void getBestRateScorePair(const MapRateToScoreType& map_rate_to_score, double& best_rate, double& best_score)
  {
    best_score = -1;
    for (MapRateToScoreType::const_iterator mit = map_rate_to_score.begin(); mit != map_rate_to_score.end(); ++mit)
    {
      if (mit->second > best_score)
      {
        best_score = mit->second;
        best_rate = mit->first;
      }
    }
  }

  PeakSpectrum extractPeakSpectrum(Size element_count, double mass_diff, double rt, double feature_hit_theoretical_mz, Int feature_hit_charge, const PeakMap& peak_map)
  {
    PeakSpectrum spec = *peak_map.RTBegin(rt - 1e-8);
    PeakSpectrum::ConstIterator begin_it = spec.MZBegin(feature_hit_theoretical_mz - 1e-8);
    PeakSpectrum::ConstIterator end_it = spec.MZEnd(feature_hit_theoretical_mz + element_count * mass_diff / feature_hit_charge + 1e-8);

    PeakSpectrum ret;
    for (; begin_it != end_it; ++begin_it)
    {
      if (begin_it->getIntensity() > 1e-8)
      {
        ret.push_back(*begin_it);
      }
    }
    return ret;
  }

  // collects intensities starting at seed_mz/_rt, if no peak is found at the expected position a 0 is added
  vector<double> extractIsotopicIntensities(Size element_count, double mass_diff, double mz_tolerance_ppm,
                                            double seed_rt, double seed_mz, double seed_charge,
                                            const PeakMap& peak_map)
  {
    vector<double> isotopic_intensities;
    for (Size k = 0; k != element_count; ++k)
    {
      double min_rt = seed_rt - 0.01; // feature rt
      double max_rt = seed_rt + 0.01;
      double mz = seed_mz + k * mass_diff / seed_charge;

      double min_mz;
      double max_mz;

      if (k <= 5)
      {
        double ppm = std::max(10.0, mz_tolerance_ppm); // restrict ppm to 10 for low intensity peaks
        min_mz = mz - mz * ppm * 1e-6;
        max_mz = mz + mz * ppm * 1e-6;
      }
      else
      {
        min_mz = mz - mz * mz_tolerance_ppm * 1e-6;
        max_mz = mz + mz * mz_tolerance_ppm * 1e-6;
      }

      double found_peak_int = 0;

      PeakMap::ConstAreaIterator aait = peak_map.areaBeginConst(min_rt, max_rt, min_mz, max_mz);

      // find 13C/15N peak in window around theoretical predicted position
      vector<double> found_peaks;
      for (; aait != peak_map.areaEndConst(); ++aait)
      {
        double peak_int = aait->getIntensity();
        if (peak_int > 1) // we found a valid 13C/15N peak
        {
          found_peaks.push_back(peak_int);
        }
      }

      found_peak_int = std::accumulate(found_peaks.begin(), found_peaks.end(), 0.0);

      // assign peak intensity to first peak in small area around theoretical predicted position (should be usually only be 1)
      isotopic_intensities.push_back(found_peak_int);
    }
    return isotopic_intensities;
  }

  void writePeakIntensities_(SVOutStream& out_stream, vector<double> isotopic_intensities, bool write_13Cpeaks)
  {
    double intensities_sum_12C = 0.0;
    // calculate 12C summed intensity
    for (Size k = 0; k != 5; ++k)
    {
      if (k >= isotopic_intensities.size())
      {
        break;
      }
      intensities_sum_12C += isotopic_intensities[k];
    }

    // determine 13C peaks and summed intensity
    double intensities_sum_13C = 0;
    for (Size u = 5; u < isotopic_intensities.size(); ++u)
    {
      intensities_sum_13C += isotopic_intensities[u];
    }

    String int_string;
    // print 12C peaks
    for (Size u = 0; u != 5; ++u)
    {
      if (u == isotopic_intensities.size())
      {
        break;
      }
      int_string += String::number(isotopic_intensities[u], 0);
      int_string += " ";
    }
    int_string += ", ";

    if (write_13Cpeaks)
    {
      // print 13C peaks
      for (Size u = 5; u < isotopic_intensities.size(); ++u)
      {
        int_string += String::number(isotopic_intensities[u], 0);
        if (u < isotopic_intensities.size() - 1)
        {
          int_string += " ";
        }
      }
      out_stream << int_string;

      double ratio = 0.0;
      if (intensities_sum_12C + intensities_sum_13C > 0.0000001)
      {
        ratio = intensities_sum_13C / (intensities_sum_12C + intensities_sum_13C);
      }
      out_stream << ratio; // << skewness(intensities_13C.begin(), intensities_13C.end());

    }
    else // bad correlation, no need to print intensities, ratio etc.
    {
      out_stream << "\t\t";
    }
  }

  // scores smaller than 0 will be padded to 0
  MapRateToScoreType normalizeToMax(const MapRateToScoreType& map_rate_to_decomposition_weight)
  {
    // extract highest weight (best score) and rate
    double best_rate, best_score;
    getBestRateScorePair(map_rate_to_decomposition_weight, best_rate, best_score);

    if (debug_level_ >= 10)
    {
      OPENMS_LOG_DEBUG << "best rate + score: " << best_rate << " " << best_score << endl;
    }

    // normalize weights to max(weights)=1
    MapRateToScoreType map_weights_norm(map_rate_to_decomposition_weight);
    for (MapRateToScoreType::iterator mit = map_weights_norm.begin(); mit != map_weights_norm.end(); ++mit)
    {
      if (best_score > 0)
      {
        mit->second /= best_score;
      }
      else
      {
        mit->second = 0;
      }
    }

    return map_weights_norm;
  }

  // Extract the mono-isotopic trace and reports the rt of the maximum intensity
  // Used to compensate for slight RT shifts (e.g. important if features of a different map are used)
  // n_scans corresponds to the number of neighboring scan rts that should be extracted
  // n_scan = 2 -> vector size = 1 + 2 + 2
  vector<double> findApexRT(const FeatureMap::iterator feature_it, double hit_rt, const PeakMap& peak_map, Size n_scans)
  {
    vector<double> seeds_rt;
    vector<Peak2D> mono_trace;

    if (!feature_it->getConvexHulls().empty())
    {
      // extract elution profile of 12C containing mass trace using a bounding box
      // first convex hull contains the monoisotopic 12C trace
      const DBoundingBox<2>& mono_bb = feature_it->getConvexHulls()[0].getBoundingBox();

      //(min_rt, max_rt, min_mz, max_mz)
      PeakMap::ConstAreaIterator ait = peak_map.areaBeginConst(mono_bb.minPosition()[0], mono_bb.maxPosition()[0], mono_bb.minPosition()[1], mono_bb.maxPosition()[1]);
      for (; ait != peak_map.areaEndConst(); ++ait)
      {
        Peak2D p2d;
        p2d.setRT(ait.getRT()); // get rt of scan
        p2d.setMZ(ait->getMZ()); // get peak 1D mz
        p2d.setIntensity(ait->getIntensity());
        mono_trace.push_back(p2d);
      }
    }

    // if there is no 12C mono trace generate a valid starting point
    if (mono_trace.empty())
    {
      Peak2D p2d;
      double next_valid_scan_rt = peak_map.RTBegin(hit_rt - 0.001)->getRT();
      p2d.setRT(next_valid_scan_rt);
      p2d.setMZ(0); // actually not needed
      p2d.setIntensity(0);
      mono_trace.push_back(p2d);
    }

    // determine trace peak with highest intensity
    double max_trace_int = -1e16;
    Size max_trace_int_idx = 0;

    for (Size j = 0; j != mono_trace.size(); ++j)
    {
      if (mono_trace[j].getIntensity() > max_trace_int)
      {
        max_trace_int = mono_trace[j].getIntensity();
        max_trace_int_idx = j;
      }
    }
    double max_trace_int_rt = mono_trace[max_trace_int_idx].getRT();
    seeds_rt.push_back(max_trace_int_rt);

    for (Size i = 1; i <= n_scans; ++i)
    {
      double rt_after = max_trace_int_rt;
      if (max_trace_int_idx < mono_trace.size() - (Int)i)
      {
        rt_after = mono_trace[max_trace_int_idx + i].getRT();
      }

      double rt_before = max_trace_int_rt;
      if (max_trace_int_idx >= i)
      {
        rt_before = mono_trace[max_trace_int_idx - i].getRT();
      }

      if (fabs(max_trace_int_rt - rt_after) < 10.0)
      {
        seeds_rt.push_back(rt_after);
      }

      if (fabs(max_trace_int_rt - rt_before) < 10.0)
      {
        seeds_rt.push_back(rt_before);
      }
    }
    //cout << "Seeds size:" << seeds_rt.size() << endl;
    return seeds_rt;
  }

  PeakSpectrum mergeSpectra(const PeakMap& to_merge)
  {
    PeakSpectrum merged;
    for (Size i = 0; i != to_merge.size(); ++i)
    {
      std::copy(to_merge[i].begin(), to_merge[i].end(), std::back_inserter(merged));
    }
    merged.sortByPosition();

    return merged;
  }

  ///> converts a vector of isotopic intensities to a peak spectrum starting at mz=mz_start with mass_diff/charge step size
  PeakSpectrum isotopicIntensitiesToSpectrum(double mz_start, double mass_diff, Int charge, vector<double> isotopic_intensities)
  {
    PeakSpectrum ps;
    for (Size i = 0; i != isotopic_intensities.size(); ++i)
    {
      Peak1D peak;
      peak.setMZ(mz_start + i * mass_diff / (double)charge);
      peak.setIntensity(isotopic_intensities[i]);
      ps.push_back(peak);
    }
    return ps;
  }

  ///> Collect decomposition coefficients in the merge window around the correlation maximum.
  ///> Final list of RIAs is constructed for the peptide.
  void extractIncorporationsAtCorrelationMaxima(SIPPeptide& sip_peptide,
                                                const IsotopePatterns& patterns,
                                                double weight_merge_window = 5.0,
                                                double min_corr_threshold = 0.5,
                                                double min_decomposition_weight = 10.0)
  {
    const MapRateToScoreType& map_rate_to_decomposition_weight = sip_peptide.decomposition_map;
    const MapRateToScoreType& map_rate_to_correlation_score = sip_peptide.correlation_map;
    vector<SIPIncorporation> sip_incorporations;
    const vector<RateScorePair>& corr_maxima = sip_peptide.correlation_maxima;

    double explained_TIC_fraction = 0;
    double TIC = 0;
    Size non_zero_decomposition_coefficients = 0;

    double max_corr_TIC = 0;

    for (Size k = 0; k < corr_maxima.size(); ++k)
    {
      const double rate = corr_maxima[k].rate;
      const double corr = corr_maxima[k].score;

      if (corr > min_corr_threshold)
      {
        SIPIncorporation sip_incorporation{};
        sip_incorporation.rate = rate;

        // sum up decomposition intensities for quantification in merge window
        double int_sum = 0;
        MapRateToScoreType::const_iterator low = map_rate_to_decomposition_weight.lower_bound(rate - weight_merge_window - 1e-4);
        MapRateToScoreType::const_iterator high = map_rate_to_decomposition_weight.lower_bound(rate + weight_merge_window + 1e-4);
        for (; low != high; ++low)
        {
          int_sum += low->second;
        }

        if (low != map_rate_to_decomposition_weight.end())
        {
          int_sum += low->second;
        }

        sip_incorporation.abundance = int_sum; // calculate abundance as sum of all decompositions
        sip_incorporation.correlation = min(corr, 1.0);

        max_corr_TIC += int_sum;

        // find closest idx (could be more efficient using binary search)
        Size closest_idx = 0;
        for (Size i = 0; i != patterns.size(); ++i)
        {
          if (fabs(patterns[i].first - rate) < fabs(patterns[closest_idx].first - rate))
          {
            closest_idx = i;
          }
        }
#ifdef DEBUG_METAPROSIP
        sip_incorporation.theoretical = isotopicIntensitiesToSpectrum(sip_peptide.mz_theo, sip_peptide.mass_diff, sip_peptide.charge, patterns[closest_idx].second);
#endif
        if (int_sum > 1e-4)
        {
          sip_incorporations.push_back(sip_incorporation);
        }
        else
        {
          if (debug_level_ > 1)
          {
            OPENMS_LOG_WARN << "warning: prevented adding of 0 abundance decomposition at rate " << rate << endl;
            OPENMS_LOG_WARN << "decomposition: " << endl;
            for (MapRateToScoreType::const_iterator it = map_rate_to_decomposition_weight.begin(); it != map_rate_to_decomposition_weight.end(); ++it)
            {
              OPENMS_LOG_WARN << it->first << " " << it->second << endl;
            }
            OPENMS_LOG_WARN << "correlation: " << endl;
            for (MapRateToScoreType::const_iterator it = map_rate_to_correlation_score.begin(); it != map_rate_to_correlation_score.end(); ++it)
            {
              OPENMS_LOG_WARN << it->first << " " << it->second << endl;
            }
          }

        }
      }
    }

    // find highest non-natural incorporation
    double highest_non_natural_abundance = 0;
    double highest_non_natural_rate = 0;
    for (vector<SIPIncorporation>::const_iterator it = sip_incorporations.begin(); it != sip_incorporations.end(); ++it)
    {
      if (it->rate < 5.0) // skip natural
      {
        continue;
      }

      if (it->abundance > highest_non_natural_abundance)
      {
        highest_non_natural_rate = it->rate;
        highest_non_natural_abundance = it->abundance;
      }
    }

    bool non_natural = false;
    if (highest_non_natural_rate > 5.0 && highest_non_natural_abundance > min_decomposition_weight)
    {
      non_natural = true;
    }

    // used for non-gaussian shape detection
    for (MapRateToScoreType::const_iterator mit = map_rate_to_decomposition_weight.begin(); mit != map_rate_to_decomposition_weight.end(); ++mit)
    {
      double decomposition_rate = mit->first;
      double decomposition_weight = mit->second;
      TIC += decomposition_weight;

      if (non_natural && decomposition_weight > 0.05 * highest_non_natural_abundance && decomposition_rate > 5.0)
      {
        ++non_zero_decomposition_coefficients;
      }
    }

    if (TIC > 1e-5)
    {
      explained_TIC_fraction = max_corr_TIC / TIC;
    }
    else
    {
      explained_TIC_fraction = 0;
    }

    // set results
    sip_peptide.incorporations = sip_incorporations;
    sip_peptide.explained_TIC_fraction = explained_TIC_fraction;
    sip_peptide.non_zero_decomposition_coefficients = non_zero_decomposition_coefficients;
  }

  ///> Collect decomposition coefficients. Starting at the largest decomposition weights merge smaller weights in the merge window.
  void extractIncorporationsAtHeighestDecompositionWeights(SIPPeptide& sip_peptide,
                                                           const IsotopePatterns& patterns,
                                                           double weight_merge_window = 5.0,
                                                           double min_corr_threshold = 0.5,
                                                           double min_low_RIA_threshold = -1,
                                                           double min_decomposition_weight = 10.0)
  {
    if (min_low_RIA_threshold < 0)
    {
      min_low_RIA_threshold = min_corr_threshold;
    }

    const MapRateToScoreType& map_rate_to_decomposition_weight = sip_peptide.decomposition_map;
    const MapRateToScoreType& map_rate_to_correlation_score = sip_peptide.correlation_map;

    double explained_TIC_fraction = 0;
    double TIC = 0;
    Size non_zero_decomposition_coefficients = 0;
    double max_corr_TIC = 0;
    vector<SIPIncorporation> sip_incorporations;

    // find decomposition weights with correlation larger than threshold (seeds)
    MapRateToScoreType::const_iterator md_it = map_rate_to_decomposition_weight.begin();
    MapRateToScoreType::const_iterator mc_it = map_rate_to_correlation_score.begin();

    set<pair<double, double> > seeds_weight_rate_pair;
    for (; md_it != map_rate_to_decomposition_weight.end(); ++md_it, ++mc_it)
    {
      if (mc_it->first < 10.0) // lowRIA region
      {
        if (mc_it->second >= min_low_RIA_threshold && md_it->second >= min_decomposition_weight)
        {
          seeds_weight_rate_pair.insert(make_pair(md_it->second, md_it->first));
        }
      }
      else // non-low RIA region
      {
        if (mc_it->second >= min_corr_threshold && md_it->second >= min_decomposition_weight)
        {
          seeds_weight_rate_pair.insert(make_pair(md_it->second, md_it->first));
          //cout << "Seeds insert: " << md_it->second << " " << md_it->first << endl;
        }
      }
    }

    // cout << "Seeds: " << seeds_weight_rate_pair.size() << endl;

    // seeds_weight_rate_pair contains the seeds ordered by their decomposition weight
    while (!seeds_weight_rate_pair.empty())
    {
      // pop last element from set
      set<pair<double, double> >::iterator last_element = --seeds_weight_rate_pair.end();
      pair<double, double> current_seed = *last_element;

      //cout << current_seed.first << " " << current_seed.second << endl;

      // find weights in window to merge, remove from seed map. maybe also remove from original map depending on whether we want to quantify the weight only 1 time
      const double rate = current_seed.second;

      SIPIncorporation sip_incorporation{};
      sip_incorporation.rate = rate;

      MapRateToScoreType::const_iterator low = map_rate_to_decomposition_weight.lower_bound(rate - weight_merge_window - 1e-4);
      MapRateToScoreType::const_iterator high = map_rate_to_decomposition_weight.lower_bound(rate + weight_merge_window + 1e-4);

      // cout << "Distance: " << std::distance(low, high) << endl;;

      MapRateToScoreType::const_iterator l1 = low;
      MapRateToScoreType::const_iterator h1 = high;

      // iterate over peaks in merge window
      for (; l1 != h1; ++l1)
      {
        // remove from seed map
        seeds_weight_rate_pair.erase(make_pair(l1->second, l1->first));
      }

      // Sum up decomposition intensities for quantification in merge window
      double int_sum = 0;
      for (; low != high; ++low)
      {
        int_sum += low->second;
      }

      if (low != map_rate_to_decomposition_weight.end())
      {
        int_sum += low->second;
      }

      sip_incorporation.abundance = int_sum;
      MapRateToScoreType::const_iterator corr_it = map_rate_to_correlation_score.lower_bound(rate - 1e-6);
      sip_incorporation.correlation = min(corr_it->second, 1.0);

      max_corr_TIC += int_sum;

      PeakSpectrum theoretical_spectrum;

      // find closest idx (could be more efficient using binary search)
      Size closest_idx = 0;
      for (Size i = 0; i != patterns.size(); ++i)
      {
        if (fabs(patterns[i].first - rate) < fabs(patterns[closest_idx].first - rate))
        {
          closest_idx = i;
        }
      }

#ifdef DEBUG_METAPROSIP
      sip_incorporation.theoretical = isotopicIntensitiesToSpectrum(sip_peptide.mz_theo, sip_peptide.mass_diff, sip_peptide.charge, patterns[closest_idx].second);
#endif

      sip_incorporations.push_back(sip_incorporation);
    }

    // find highest non-natural incorporation
    double highest_non_natural_abundance = 0;
    double highest_non_natural_rate = 0;
    for (vector<SIPIncorporation>::const_iterator it = sip_incorporations.begin(); it != sip_incorporations.end(); ++it)
    {
      if (it->rate < 5.0) // skip natural
      {
        continue;
      }

      if (it->abundance > highest_non_natural_abundance)
      {
        highest_non_natural_rate = it->rate;
        highest_non_natural_abundance = it->abundance;
      }
    }

    bool non_natural = false;
    if (highest_non_natural_rate > 5.0)
    {
      non_natural = true;
    }

    // used for non-gaussian shape detection
    for (MapRateToScoreType::const_iterator mit = map_rate_to_decomposition_weight.begin(); mit != map_rate_to_decomposition_weight.end(); ++mit)
    {
      double decomposition_rate = mit->first;
      double decomposition_weight = mit->second;
      TIC += decomposition_weight;

      if (non_natural && decomposition_weight > 0.05 * highest_non_natural_abundance && decomposition_rate > 5.0)
      {
        ++non_zero_decomposition_coefficients;
      }
    }

    if (TIC > 1e-5)
    {
      explained_TIC_fraction = max_corr_TIC / TIC;
    }
    else
    {
      explained_TIC_fraction = 0;
    }

    // set results
    std::sort(sip_incorporations.begin(), sip_incorporations.end(), RIALess());
    sip_peptide.incorporations = sip_incorporations;
    sip_peptide.explained_TIC_fraction = explained_TIC_fraction;
    sip_peptide.non_zero_decomposition_coefficients = non_zero_decomposition_coefficients;
  }

  ///> calculate the global labeling ration based on all but the first 4 peaks
  double calculateGlobalLR(const vector<double>& isotopic_intensities)
  {
    if (isotopic_intensities.size() < 5)
    {
      return 0.0;
    }

    double sum = accumulate(isotopic_intensities.begin(), isotopic_intensities.end(), 0.0);
    double sum_incorporated = accumulate(isotopic_intensities.begin() + 4, isotopic_intensities.end(), 0.0);

    if (sum < 1e-4)
    {
      return 0.0;
    }

    return sum_incorporated / sum;
  }

  ExitCodes main_(int, const char**) override
  {
    String file_extension_ = getStringOption_("plot_extension");
    Int debug_level = getIntOption_("debug");
    String in_mzml = getStringOption_("in_mzML");
    String in_features = getStringOption_("in_featureXML");
    double mz_tolerance_ppm_ = getDoubleOption_("mz_tolerance_ppm");
    double rt_tolerance_s = getDoubleOption_("rt_tolerance_s");

    double weight_merge_window_ = getDoubleOption_("weight_merge_window");
    double intensity_threshold_ = getDoubleOption_("intensity_threshold");
    double decomposition_threshold = getDoubleOption_("decomposition_threshold");

    Size min_consecutive_isotopes = (Size)getIntOption_("min_consecutive_isotopes");

    String qc_output_directory = getStringOption_("qc_output_directory");

    Size n_heatmap_bins = getIntOption_("heatmap_bins");
    double score_plot_y_axis_min = getDoubleOption_("score_plot_yaxis_min");

    String tmp_path = File::getTempDirectory();
    tmp_path.substitute('\\', '/');

    // Do we want to create a qc report?
    if (!qc_output_directory.empty())
    {
      QString executable = getStringOption_("r_executable").toQString();
      // convert path to absolute path
      QDir qc_dir(qc_output_directory.toQString());
      qc_output_directory = String(qc_dir.absolutePath());

      // trying to create qc_output_directory if not present
      if (!qc_dir.exists())
      {
        qc_dir.mkpath(qc_output_directory.toQString());
      }
      // check if R and dependencies are installed
      StringList package_names;
      package_names.push_back("gplots");

      bool R_is_working = RIntegration::checkRDependencies(tmp_path, package_names, executable);
      if (!R_is_working)
      {
        OPENMS_LOG_INFO << "There was a problem detecting one of the required R libraries." << endl;
        return EXTERNAL_PROGRAM_ERROR;
      }
    }

    String out_csv = getStringOption_("out_csv");
    ofstream out_csv_stream(out_csv.c_str());
    out_csv_stream << fixed << setprecision(4);

    String out_peptide_centric_csv = getStringOption_("out_peptide_centric_csv");
    ofstream out_peptide_csv_stream(out_peptide_centric_csv.c_str());
    out_peptide_csv_stream << fixed << setprecision(4);

    String labeling_element = getStringOption_("labeling_element");

    //bool plot_merged = getFlag_("plot_merged");
    bool report_natural_peptides = getFlag_("report_natural_peptides");
    bool use_unassigned_ids = getFlag_("use_unassigned_ids");
    bool use_averagine_ids = getFlag_("use_averagine_ids");

    //String debug_patterns_name = getStringOption_("debug_patterns_name");

    double correlation_threshold = getDoubleOption_("correlation_threshold");

    double xic_threshold = getDoubleOption_("xic_threshold");

    double min_correlation_distance_to_averagine = getDoubleOption_("min_correlation_distance_to_averagine");

    bool cluster_flag = getFlag_("cluster");

    // read descriptions from FASTA and create map for fast annotation
    String in_fasta = getStringOption_("in_fasta");
    vector<FASTAFile::FASTAEntry> fasta_entries;
    FASTAFile fasta_file;
    fasta_file.setLogType(log_type_);
    fasta_file.load(in_fasta, fasta_entries);
    map<String, String> proteinid_to_description;
    for (vector<FASTAFile::FASTAEntry>::const_iterator it = fasta_entries.begin(); it != fasta_entries.end(); ++it)
    {
      if (!it->identifier.empty() && !it->description.empty())
      {
        String s = it->identifier;
        proteinid_to_description[s.trim().toUpper()] = it->description;
      }
    }

    OPENMS_LOG_INFO << "loading feature map..." << endl;
    FeatureMap feature_map;
    FileHandler().loadFeatures(in_features, feature_map, {FileTypes::FEATUREXML});

    // annotate as features found using feature finding (to distinguish them from averagine features oder id based features ... see below)
    for (FeatureMap::iterator feature_it = feature_map.begin(); feature_it != feature_map.end(); ++feature_it)
    {
      feature_it->setMetaValue("feature_type", FEATURE_STRING);
    }

    // if also unassigned ids are used create a pseudo feature
    if (use_unassigned_ids)
    {
      const vector<PeptideIdentification> unassigned_ids = feature_map.getUnassignedPeptideIdentifications();
      Size unassigned_id_features = 0;
      for (vector<PeptideIdentification>::const_iterator it = unassigned_ids.begin(); it != unassigned_ids.end(); ++it)
      {
        vector<PeptideHit> hits = it->getHits();
        if (!hits.empty())
        {
          Feature f;
          f.setMetaValue("feature_type", UNASSIGNED_ID_STRING);
          f.setRT(it->getRT());
          // take sequence of first hit to calculate ground truth mz
          Int charge = hits[0].getCharge();
          if (charge == 0)
          {
            continue;
          }
          double mz =  hits[0].getSequence().getMZ(charge);
          f.setMZ(mz);
          // add id to pseudo feature
          vector<PeptideIdentification> id;
          id.push_back(*it);
          f.setPeptideIdentifications(id);
          feature_map.push_back(f);
          unassigned_id_features++;
        }
      }
      feature_map.updateRanges();
      OPENMS_LOG_INFO << "Evaluating " << unassigned_id_features << " unassigned identifications." << endl;
    }

    // determine all spectra that have not been identified and assign an averagine peptide to it
    if (use_averagine_ids)
    {
      // load only MS2 spectra with precursor information
      PeakMap peak_map;
      FileHandler mh;
      std::vector<Int> ms_level(1, 2);
      mh.getOptions().setMSLevels(ms_level);
      mh.loadExperiment(in_mzml, peak_map, {FileTypes::MZML});
      peak_map.sortSpectra();
      peak_map.updateRanges();

      // extract rt and mz of all identified precursors and store them in blacklist
      vector<Peak2D> blacklisted_precursors;
      // in features
      for (FeatureMap::iterator feature_it = feature_map.begin(); feature_it != feature_map.end(); ++feature_it) // for each peptide feature
      {
        const vector<PeptideIdentification>& f_ids = feature_it->getPeptideIdentifications();
        for (vector<PeptideIdentification>::const_iterator id_it = f_ids.begin(); id_it != f_ids.end(); ++id_it)
        {
          if (!id_it->getHits().empty())
          {
            // Feature with id found so we don't need to generate averagine id. Find MS2 in experiment and blacklist it.
            Peak2D p;
            p.setRT(id_it->getRT());
            p.setMZ(id_it->getMZ());
            blacklisted_precursors.push_back(p);
          }
        }
      }

      // and in unassigned ids
      const vector<PeptideIdentification> unassigned_ids = feature_map.getUnassignedPeptideIdentifications();
      for (vector<PeptideIdentification>::const_iterator it = unassigned_ids.begin(); it != unassigned_ids.end(); ++it)
      {
        const vector<PeptideHit> hits = it->getHits();
        if (!hits.empty())
        {
          Peak2D p;
          p.setRT(it->getRT());
          p.setMZ(it->getMZ());
          blacklisted_precursors.push_back(p);
        }
      }

      // find index of all precursors that have been blacklisted
      vector<Size> blacklist_idx;
      for (vector<Peak2D>::const_iterator it = blacklisted_precursors.begin(); it != blacklisted_precursors.end(); ++it)
      {
        PeakMap::const_iterator map_rt_begin = peak_map.RTBegin(-std::numeric_limits<double>::max());
        PeakMap::const_iterator rt_begin = peak_map.RTBegin(it->getRT() - 1e-5);
        Size index = std::distance(map_rt_begin, rt_begin);
        //cout << "Blacklist Index: " << index << endl;
        blacklist_idx.push_back(index);
      }

      for (Size i = 0; i != peak_map.size(); ++i)
      {
        // precursor not blacklisted?
        if (find(blacklist_idx.begin(), blacklist_idx.end(), i) == blacklist_idx.end() && !peak_map[i].getPrecursors().empty())
        {
          // store feature with id generated from averagine peptide (pseudo id)
          Feature f;

          double precursor_mz = peak_map[i].getPrecursors()[0].getMZ();
          int precursor_charge = peak_map[i].getPrecursors()[0].getCharge();
          //double precursor_mass = (double)precursor_charge * precursor_mz - (double)precursor_charge * Constants::PROTON_MASS_U;

          // add averagine id to pseudo feature
          PeptideHit pseudo_hit;

          // set peptide with lowest deviation from averagine
          pseudo_hit.setSequence(AASequence()); // set empty sequence
          pseudo_hit.setCharge(precursor_charge);
          PeptideIdentification pseudo_id;
          vector<PeptideHit> pseudo_hits;
          pseudo_hits.push_back(pseudo_hit);
          pseudo_id.setHits(pseudo_hits);
          vector<PeptideIdentification> id;
          id.push_back(pseudo_id);
          f.setPeptideIdentifications(id);
          f.setRT(peak_map[i].getRT());
          f.setMZ(precursor_mz);
          f.setMetaValue("feature_type", UNIDENTIFIED_STRING);
          feature_map.push_back(f);
        }
      }
      feature_map.updateRanges();
    }

    OPENMS_LOG_INFO << "loading experiment..." << endl;
    PeakMap peak_map;
    FileHandler mh;
    std::vector<Int> ms_level(1, 1);
    mh.getOptions().setMSLevels(ms_level);
    mh.loadExperiment(in_mzml, peak_map, {FileTypes::MZML});
    peak_map.updateRanges();
    ThresholdMower tm;
    Param tm_parameters;
    tm_parameters.setValue("threshold", intensity_threshold_);
    tm.setParameters(tm_parameters);
    tm.filterPeakMap(peak_map);
    peak_map.sortSpectra();

    // used to generate plots
    vector<String> titles;
    vector<MapRateToScoreType> weight_maps;
    vector<MapRateToScoreType> normalized_weight_maps;
    vector<MapRateToScoreType> correlation_maps;

    String file_suffix = "_" + String(QFileInfo(in_mzml.toQString()).baseName()) + "_" + String::random(4);

    vector<SIPPeptide> sip_peptides;

    Size nPSMs = 0; ///< number of PSMs. If 0 IDMapper has not been called.
    Size spectrum_with_no_isotopic_peaks(0);
    Size spectrum_with_isotopic_peaks(0);

    for (FeatureMap::iterator feature_it = feature_map.begin(); feature_it != feature_map.end(); ++feature_it) // for each peptide feature
    {
      const double feature_hit_center_rt = feature_it->getRT();

      // check if out of experiment bounds
      if (feature_hit_center_rt > peak_map.getMaxRT() || feature_hit_center_rt < peak_map.getMinRT())
      {
        continue;
      }

      // Extract 1 or more MS/MS with identifications assigned to the feature by IDMapper
      vector<PeptideIdentification> pep_ids = feature_it->getPeptideIdentifications();

      nPSMs += pep_ids.size();

      // Skip features without peptide identifications
      if (pep_ids.empty())
      {
        continue;
      }

      // add best scoring PeptideHit of all PeptideIdentifications mapping to the current feature to tmp_pepid
      PeptideIdentification tmp_pepid;
      tmp_pepid.setHigherScoreBetter(pep_ids[0].isHigherScoreBetter());
      for (Size i = 0; i != pep_ids.size(); ++i)
      {
        pep_ids[i].assignRanks();
        const vector<PeptideHit>& hits = pep_ids[i].getHits();
        if (!hits.empty())
        {
          tmp_pepid.insertHit(hits[0]);
        }
        else
        {
          OPENMS_LOG_WARN << "Empty peptide hit encountered on feature. Ignoring." << endl;
        }
      }

      tmp_pepid.assignRanks();

      SIPPeptide sip_peptide;
      sip_peptide.feature_type = feature_it->getMetaValue("feature_type"); // used to annotate feature type in reporting

      // retrieve identification information
      const PeptideHit& feature_hit = tmp_pepid.getHits()[0];
      const double feature_hit_score = feature_hit.getScore();
      const double feature_hit_center_mz = feature_it->getMZ();
      const Int feature_hit_charge = feature_hit.getCharge();

      String feature_hit_seq = "";
      double feature_hit_theoretical_mz = 0;
      AASequence feature_hit_aaseq;
      // set theoretical mz of peptide hit to:
      //   mz of sequence if we have a sequence identified
      // otherwise:
      //   mz of precursor (stored in feature mz) if no sequence identified
      if (sip_peptide.feature_type == FEATURE_STRING || sip_peptide.feature_type == UNASSIGNED_ID_STRING)
      {
        feature_hit_aaseq = feature_hit.getSequence();
        feature_hit_seq = feature_hit_aaseq.toString();
        feature_hit_theoretical_mz = feature_hit_aaseq.getMZ(feature_hit.getCharge());
      }
      else if (sip_peptide.feature_type == UNIDENTIFIED_STRING)
      {
        feature_hit_aaseq = AASequence();
        feature_hit_seq = String("");
        feature_hit_theoretical_mz = feature_hit_center_mz;
      }

      if (debug_level_ > 1)
      {
        OPENMS_LOG_DEBUG << "Feature type: (" << sip_peptide.feature_type << ") Seq.: " << feature_hit_seq << " m/z: " << feature_hit_theoretical_mz << endl;
      }

      const set<String> protein_accessions = feature_hit.extractProteinAccessionsSet();
      sip_peptide.accessions = vector<String>(protein_accessions.begin(), protein_accessions.end());
      sip_peptide.sequence = feature_hit_aaseq;
      sip_peptide.mz_theo = feature_hit_theoretical_mz;
      sip_peptide.mass_theo = feature_hit_theoretical_mz * feature_hit_charge - feature_hit_charge * Constants::PROTON_MASS_U;
      sip_peptide.charge = feature_hit_charge;
      sip_peptide.score = feature_hit_score;
      sip_peptide.feature_rt = feature_hit_center_rt;
      sip_peptide.feature_mz = feature_hit_center_mz;
      sip_peptide.unique = sip_peptide.accessions.size() == 1;

      // determine retention time of scans next to the central scan
      vector<double> seeds_rt = findApexRT(feature_it, feature_hit_center_rt, peak_map, 2); // 1 scan at maximum, 2+2 above and below
      double max_trace_int_rt = seeds_rt[0];

      // determine maximum number of peaks and mass difference
      EmpiricalFormula e = feature_hit_aaseq.getFormula();

      // assign mass difference between labeling element isotopes
      if (labeling_element == "C")
      {
        sip_peptide.mass_diff = 1.003354837810;
      }
      else if (labeling_element == "N")
      {
        sip_peptide.mass_diff = 0.9970349;
      }
      else if (labeling_element == "H")
      {
        sip_peptide.mass_diff = 1.00627675;
      }
      else if (labeling_element == "O")
      {
        // 18O-16O distance is approx. 2.0042548 Dalton but natural isotopic pattern is dominated by 13C-12C distance (approx. 1.0033548)
        // After the convolution of the O-isotope distribution with the natural one we get multiple copies of the O-distribution (with 2 Da spaces)
        // shifted by 13C-12C distances. Choosing (18O-16O) / 2 as expected mass trace distance should therefor collect all of them.
        sip_peptide.mass_diff = 2.0042548 / 2.0;
      }

      Size element_count(0);
      Size isotopic_trace_count(0);
      if (sip_peptide.feature_type == FEATURE_STRING || sip_peptide.feature_type == UNASSIGNED_ID_STRING)
      {
        element_count = MetaProSIPDecomposition::getNumberOfLabelingElements(labeling_element, feature_hit_aaseq);
      }
      else // if (sip_peptide.feature_type == UNIDENTIFIED_STRING)
      {
        // calculate number of expected labeling elements using averagine model C:4.9384 H:7.7583 N:1.3577 O:1.4773 S:0.0417 divided by average weight 111.1254
        if (labeling_element == "C")
        {
          element_count = static_cast<Size>(sip_peptide.mass_theo * 0.0444398894906044);
        }
        else if (labeling_element == "N")
        {
          element_count = static_cast<Size>(sip_peptide.mass_theo * 0.0122177302837372);
        }
        else if (labeling_element == "H")
        {
          element_count = static_cast<Size>(sip_peptide.mass_theo * 0.06981572169);
        }
        else if (labeling_element == "O")
        {
          element_count = static_cast<Size>(sip_peptide.mass_theo * 0.01329399039);
        }
      }

      isotopic_trace_count = labeling_element != "O" ? element_count : element_count * 2;

      // collect 13C / 15N peaks
      if (debug_level_ >= 10)
      {
        OPENMS_LOG_DEBUG << "Extract XICs" << endl;
      }

      vector<double> isotopic_intensities = MetaProSIPXICExtraction::extractXICsOfIsotopeTraces(isotopic_trace_count + ADDITIONAL_ISOTOPES, sip_peptide.mass_diff, mz_tolerance_ppm_, rt_tolerance_s, max_trace_int_rt, feature_hit_theoretical_mz, feature_hit_charge, peak_map, xic_threshold);

      // set intensity to zero if not enough neighboring isotopic peaks are present
      for (Size i = 0; i != isotopic_intensities.size(); ++i)
      {
        if (isotopic_intensities[i] < 1e-4) continue;
        Size consecutive_isotopes = 0;
        Size j = i;

        while (j != std::numeric_limits<Size>::max()) // unsigned type wrap-around is well defined
        {
          if (isotopic_intensities[j] <= 1e-4) break;
          ++consecutive_isotopes;
          --j;
        }
        j = i + 1;

        while (j < isotopic_intensities.size())
        {
          if (isotopic_intensities[j] <= 1e-4) break;
          ++consecutive_isotopes;
          ++j;
        }

        if (consecutive_isotopes < min_consecutive_isotopes)
        {
          isotopic_intensities[i] = 0;
        }
      }

      double TIC = accumulate(isotopic_intensities.begin(), isotopic_intensities.end(), 0.0);

      // collect 13C / 15N peaks
      if (debug_level_ >= 10)
      {
        OPENMS_LOG_DEBUG << "TIC of XICs: " << TIC << endl;
        for (Size i = 0; i != isotopic_intensities.size(); ++i)
        {
          cout << isotopic_intensities[i] << endl;
        }
      }

      // no Peaks collected
      if (TIC < 1e-4)
      {
        ++spectrum_with_no_isotopic_peaks;
        if (debug_level > 0)
        {
          OPENMS_LOG_INFO << "no isotopic peaks in spectrum" << endl;
        }
        continue;
      }
      else
      {
        ++spectrum_with_isotopic_peaks;
      }

      // store accumulated intensities at theoretical positions
      sip_peptide.accumulated = isotopicIntensitiesToSpectrum(feature_hit_theoretical_mz, sip_peptide.mass_diff, feature_hit_charge, isotopic_intensities);

      sip_peptide.global_LR = calculateGlobalLR(isotopic_intensities);

      Size non_zero_isotopic_intensities(0);
      for (Size i = 0; i != isotopic_intensities.size(); ++i)
      {
        if (isotopic_intensities[i] > 0.1)
        {
          ++non_zero_isotopic_intensities;
        }
      }

      if (debug_level > 0)
      {
        cout << "Isotopic intensities found / total: " << non_zero_isotopic_intensities << "/" << isotopic_intensities.size() << endl;
      }

      OPENMS_LOG_INFO << feature_hit.getSequence().toString() << "\trt: " << max_trace_int_rt << endl;

      // correlation filtering
      MapRateToScoreType map_rate_to_correlation_score;

      IsotopePatterns patterns;

      // calculate isotopic patterns for the given sequence, incoroporation interval/steps
      if (sip_peptide.feature_type == FEATURE_STRING || sip_peptide.feature_type == UNASSIGNED_ID_STRING)
      {
       if (labeling_element == "N")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor15NRange(AASequence::fromString(feature_hit_seq));
       }
       else if (labeling_element == "C")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor13CRange(AASequence::fromString(feature_hit_seq));
       }
       else if (labeling_element == "H")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor2HRange(AASequence::fromString(feature_hit_seq));
       }
       else if (labeling_element == "O")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor18ORange(AASequence::fromString(feature_hit_seq));
       }
      }
      else if (sip_peptide.feature_type == UNIDENTIFIED_STRING)
      {
       if (labeling_element == "N")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor15NRangeOfAveraginePeptide(sip_peptide.mass_theo);
       }
       else if (labeling_element == "C")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor13CRangeOfAveraginePeptide(sip_peptide.mass_theo);
       }
       else if (labeling_element == "H")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor2HRangeOfAveraginePeptide(sip_peptide.mass_theo);
       }
       else if (labeling_element == "O")
       {
         patterns = MetaProSIPDecomposition::calculateIsotopePatternsFor18ORangeOfAveraginePeptide(sip_peptide.mass_theo);
       }
      }

      // store theoretical patterns for visualization
      sip_peptide.patterns = patterns;
      for (IsotopePatterns::const_iterator pit = sip_peptide.patterns.begin(); pit != sip_peptide.patterns.end(); ++pit)
      {
        PeakSpectrum p = isotopicIntensitiesToSpectrum(feature_hit_theoretical_mz, sip_peptide.mass_diff, feature_hit_charge, pit->second);
        p.setMetaValue("rate", (double)pit->first);
        p.setMSLevel(2);
#ifdef DEBUG_METAPROSIP
        sip_peptide.pattern_spectra.push_back(p);
#endif
      }

      // calculate decomposition into isotopic patterns
      MapRateToScoreType map_rate_to_decomposition_weight;
      MetaProSIPDecomposition::calculateDecompositionWeightsIsotopicPatterns(isotopic_trace_count, isotopic_intensities, patterns, map_rate_to_decomposition_weight, sip_peptide);

      // set first intensity to zero and remove first 2 possible RIAs (0% and e.g. 1.07% for carbon)
      MapRateToScoreType tmp_map_rate_to_correlation_score;
      if (getFlag_("filter_monoisotopic"))
      {
        // calculate correlation of natural RIAs (for later reporting) before we subtract the intensities. This is somewhat redundant but no speed bottleneck.
        calculateCorrelation(isotopic_trace_count, isotopic_intensities, patterns, tmp_map_rate_to_correlation_score, labeling_element, sip_peptide.mass_theo, -1.0);
        for (Size i = 0; i != sip_peptide.reconstruction_monoistopic.size(); ++i)
        {
          if (i == 0)
          {
            isotopic_intensities[0] = 0;
          }

          isotopic_intensities[i] -= sip_peptide.reconstruction_monoistopic[i];
          if (isotopic_intensities[i] < 0)
          {
            isotopic_intensities[i] = 0;
          }
        }
      }

      sip_peptide.decomposition_map = map_rate_to_decomposition_weight;

      // calculate Pearson correlation coefficients
      calculateCorrelation(isotopic_trace_count, isotopic_intensities, patterns, map_rate_to_correlation_score, labeling_element, sip_peptide.mass_theo, min_correlation_distance_to_averagine);

      // restore original correlation of natural RIAs (take maximum of observed correlations)
      if (getFlag_("filter_monoisotopic"))
      {
        MapRateToScoreType::iterator dc_it = map_rate_to_correlation_score.begin();
        MapRateToScoreType::const_iterator tmp_dc_it = tmp_map_rate_to_correlation_score.begin();
        dc_it->second = max(tmp_dc_it->second, dc_it->second);
        ++dc_it;
        ++tmp_dc_it;
        dc_it->second = max(tmp_dc_it->second, dc_it->second);
      }

      sip_peptide.correlation_map = map_rate_to_correlation_score;

      // determine maximum correlations
      sip_peptide.correlation_maxima = MetaProSIPInterpolation::getHighPoints(correlation_threshold, map_rate_to_correlation_score);

      // FOR REPORTING: store incorporation information like e.g. theoretical spectrum for best correlations
      if (getStringOption_("collect_method") == "correlation_maximum")
      {
        extractIncorporationsAtCorrelationMaxima(sip_peptide, patterns, weight_merge_window_, correlation_threshold);
      }
      else if (getStringOption_("collect_method") == "decomposition_maximum")
      {
        extractIncorporationsAtHeighestDecompositionWeights(sip_peptide, patterns, weight_merge_window_, correlation_threshold, getDoubleOption_("lowRIA_correlation_threshold"));
      }

      // store sip peptide
      if (!sip_peptide.incorporations.empty() && sip_peptide.RR > decomposition_threshold)
      {
        if (debug_level > 0)
        {
          OPENMS_LOG_INFO << "SIP peptides: " << sip_peptide.incorporations.size() << endl;
        }
        sip_peptides.push_back(sip_peptide);
      }

      MapRateToScoreType map_rate_to_normalized_weight = normalizeToMax(map_rate_to_decomposition_weight);

      // store for plotting
      titles.push_back(feature_hit_seq + " " + String(feature_hit_center_rt));
      weight_maps.push_back(map_rate_to_decomposition_weight);
      normalized_weight_maps.push_back(map_rate_to_normalized_weight);
      correlation_maps.push_back(map_rate_to_correlation_score);
    }

    OPENMS_LOG_INFO << "Spectra with / without isotopic peaks " << spectrum_with_isotopic_peaks << "/" << spectrum_with_no_isotopic_peaks << endl;

    if (nPSMs == 0)
    {
      OPENMS_LOG_ERROR << "No assigned identifications found in featureXML. Did you forget to run IDMapper?" << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    if (sip_peptides.empty())
    {
      OPENMS_LOG_ERROR << "No peptides passing the incorporation threshold found." << endl;
      return INCOMPATIBLE_INPUT_DATA;
    }

    // copy meta information
    PeakMap debug_exp = peak_map;
    debug_exp.clear(false);

    vector<vector<SIPPeptide> > sippeptide_clusters; // vector of cluster

    if (cluster_flag)
    {
      if (debug_level > 0)
      {
        OPENMS_LOG_INFO << "Determine cluster center of RIAs: " << endl;
      }
      vector<double> cluster_center(MetaProSIPClustering::getRIAClusterCenter(sip_peptides));
      if (debug_level > 0)
      {
        OPENMS_LOG_INFO << "Assigning peptides to cluster: " << endl;
      }
      sippeptide_clusters = MetaProSIPClustering::clusterSIPPeptides(cluster_center, sip_peptides);

      // remove cluster with no assigned SIP peptide (spurious highpoints giving rise to cluster may happen because of small bumps caused by interpolation)
      vector<vector<SIPPeptide> >::iterator scit = sippeptide_clusters.begin();
      vector<double>::iterator ccit = cluster_center.begin();
      while (scit != sippeptide_clusters.end() && ccit != cluster_center.end())
      {
        if (scit->empty())
        {
          scit = sippeptide_clusters.erase(scit); // remove cluster of SIP peptides
          ccit = cluster_center.erase(ccit); // remove cluster center
        }
        else
        {
          ++scit;
          ++ccit;
        }
      }

      if (debug_level > 0)
      {
        for (Size i = 0; i != sippeptide_clusters.size(); ++i)
        {
          OPENMS_LOG_INFO << "Cluster: " << (i + 1) << " contains " << sippeptide_clusters[i].size() << " peptides." << endl;
        }
      }
    }
    else // data hasn't been clustered so just add all SIP peptides as cluster zero
    {
      sippeptide_clusters.push_back(sip_peptides);
    }

    // create group/cluster centric report
    if (!out_csv.empty())
    {
      OPENMS_LOG_INFO << "Create CSV report." << endl;
      MetaProSIPReporting::createCSVReport(sippeptide_clusters, out_csv_stream, proteinid_to_description);
    }

    // create peptide centric report
    if (!out_peptide_centric_csv.empty())
    {
      OPENMS_LOG_INFO << "Creating peptide centric report: " << out_peptide_centric_csv << std::endl;

      if (getFlag_("test"))
      {
        MetaProSIPReporting::createPeptideCentricCSVReport("test_mode_enabled.mzML", file_extension_, sippeptide_clusters, out_peptide_csv_stream, proteinid_to_description, qc_output_directory, file_suffix, report_natural_peptides);
      }
      else
      {
        MetaProSIPReporting::createPeptideCentricCSVReport(in_mzml, file_extension_, sippeptide_clusters, out_peptide_csv_stream, proteinid_to_description, qc_output_directory, file_suffix, report_natural_peptides);
      }
    }

    // quality report
    if (!qc_output_directory.empty())
    {
      QString executable = getStringOption_("r_executable").toQString();
      // TODO plot merged is now passed as false
      MetaProSIPReporting::createQualityReport(tmp_path, qc_output_directory, file_suffix, file_extension_, sippeptide_clusters, n_heatmap_bins, score_plot_y_axis_min, report_natural_peptides, executable);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  MetaProSIP tool;
  return tool.main(argc, argv);
}

///@endcond
