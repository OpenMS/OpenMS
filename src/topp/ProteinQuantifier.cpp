// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/ID/PeptideProteinResolution.h>
#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/PROCESSING/ID/IDFilter.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <cmath>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_ProteinQuantifier ProteinQuantifier

@brief Compute peptide and protein abundances from annotated feature/consensus maps or from identification results.

<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; ProteinQuantifier &rarr;</td>
            <th ALIGN = "center"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> external tools @n e.g. for statistical analysis</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
        </tr>
    </table>
</CENTER>

Reference:\n
Weisser <em>et al.</em>: <a href="https://doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

<B>Input: featureXML or consensusXML</B>

Quantification is based on the intensity values of the features in the input files. Feature intensities are first accumulated to peptide abundances, according to the peptide identifications annotated to the features/feature groups. Then, abundances of the peptides of a protein are aggregated to compute the protein abundance.

The peptide-to-protein step uses the (e.g. 3) most abundant proteotypic peptides per protein to compute the protein abundances. This is a general version of the "top 3 approach" (but only for relative quantification) described in:\n
Silva <em>et al.</em>: Absolute quantification of proteins by LCMS<sup>E</sup>: a virtue of parallel MS acquisition (Mol. Cell. Proteomics, 2006, PMID: 16219938).

Only features/feature groups with unambiguous peptide annotation are used for peptide quantification. It is possible to resolve ambiguities before applying ProteinQuantifier using one of several equivalent mechanisms in OpenMS: @ref TOPP_IDConflictResolver, @ref TOPP_ConsensusID (algorithm @p best), or @ref TOPP_FileFilter (option @p id:keep_best_score_id).

Similarly, only proteotypic peptides (i.e. those matching to exactly one protein) are used for protein quantification <em>by default</em>. Peptide/protein IDs from multiple identification runs can be handled, but will not be differentiated (i.e. protein accessions for a peptide will be accumulated over all identification runs). See section "Optional input: Protein inference/grouping results" below for exceptions to this.

Peptides with the same sequence, but with different modifications are quantified separately on the peptide level, but treated as one peptide for the protein quantification (i.e. the contributions of differently-modified variants of the same peptide are accumulated).

<B>Input: idXML</B>

Quantification based on identification results uses spectral counting, i.e. the abundance of each peptide is the number of times that peptide was identified from an MS2 spectrum (considering only the best hit per spectrum). Different identification runs in the input are treated as different samples; this makes it possible to quantify several related samples at once by merging the corresponding idXML files with @ref TOPP_IDMerger. Depending on the presence of multiple runs, output format and applicable parameters are the same as for featureXML and consensusXML, respectively.

The notes above regarding quantification on the protein level and the treatment of modifications also apply to idXML input. In particular, this means that the settings @p top 0 and @p aggregate @p sum should be used to get the "classical" spectral counting quantification on the protein level (where all identifications of all peptides of a protein are summed up).

<B>Optional input: Protein inference/grouping results</B>

By default only proteotypic peptides (i.e. those matching to exactly one protein) are used for protein quantification. However, this limitation can be overcome: Protein inference results for the whole sample set can be supplied with the @p protein_groups option (or included in a featureXML input). In that case, the peptide-to-protein references from that file are used (rather than those from @p in), and groups of indistinguishable proteins will be quantified. Each reported protein quantity then refers to the total for the respective group.

In order for everything to work correctly, it is important that the protein inference results come from the same identifications that were used to annotate the quantitative data. We suggest to use the OpenMS tool ProteinInference @ref TOPP_ProteinInference. 

More information below the parameter specification.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_ProteinQuantifier.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_ProteinQuantifier.html

<B>Output format</B>

The output files produced by this tool have a table format, with columns as described below:

<b>Protein output</b> (one protein/set of indistinguishable proteins per line):
- @b protein: Protein accession(s) (as in the annotations in the input file; separated by "/" if more than one).
- @b n_proteins: Number of indistinguishable proteins quantified (usually "1").
- @b protein_score: Protein score, e.g. ProteinProphet probability (if available).
- @b n_peptides: Number of proteotypic peptides observed for this protein (or group of indistinguishable proteins) across all samples. Note that not necessarily all of these peptides contribute to the protein abundance (depending on parameter @p top).
- @b abundance: Computed protein abundance. For consensusXML input, there will be one column  per sample ("abundance_1", "abundance_2", etc.).

<b>Peptide output</b> (one peptide or - if @p best_charge_and_fraction is set - one charge state and fraction of a peptide per line):
- @b peptide: Peptide sequence. Only peptides that occur in unambiguous annotations of features are reported.
- @b protein: Protein accession(s) for the peptide (separated by "/" if more than one).
- @b n_proteins: Number of proteins this peptide maps to. (Same as the number of accessions in the previous column.)
- @b charge: Charge state quantified in this line. "0" (for "all charges") unless @p best_charge_and_fraction was set.
- @b abundance: Computed abundance for this peptide. If the charge in the preceding column is 0, this is the total abundance of the peptide over all charge states; otherwise, it is only the abundance observed for the indicated charge (in this case, there may be more than one line for the peptide sequence). Again, for consensusXML input, there will be one column  per sample ("abundance_1", "abundance_2", etc.). Also for consensusXML, the reported values are already normalized if @p consensus:normalize was set.

<B>Protein quantification examples</B>

While quantification on the peptide level is fairly straight-forward, a number of options influence quantification on the protein level - especially for consensusXML input. The three parameters @p top:N, @p top:include_all and @p consensus:fix_peptides determine which peptides are used to quantify proteins in different samples.

As an example, consider a protein with four proteotypic peptides. Each peptide is detected in a subset of three samples, as indicated in the table below. The peptides are ranked by abundance (1: highest, 4: lowest; assuming for simplicity that the order is the same in all samples).

<CENTER>
    <table>
        <tr>
            <td></td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 1 </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 2 </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 3 </td>
        </tr>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 1 </td>
            <td ALIGN="center"> X </td>
            <td></td>
            <td ALIGN="center"> X </td>
        </tr>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 2 </td>
            <td ALIGN="center"> X </td>
            <td ALIGN="center"> X </td>
            <td></td>
        </tr>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 3 </td>
            <td ALIGN="center"> X </td>
            <td ALIGN="center"> X </td>
            <td ALIGN="center"> X </td>
        </tr>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 4 </td>
            <td ALIGN="center"> X </td>
            <td ALIGN="center"> X </td>
            <td></td>
        </tr>
    </table>
</CENTER>

Different parameter combinations lead to different quantification scenarios, as shown here:

<CENTER>
    <table>
        <tr>
          <td ALIGN="center" BGCOLOR="#EBEBEB" COLSPAN=3> @b parameters \n "*": no effect in this case </td>
          <td ALIGN="center" BGCOLOR="#EBEBEB" COLSPAN=3> <b>peptides used for quantification</b> \n "(...)": not quantified here because ... </td>
          <td ALIGN="center" VALIGN="middle" BGCOLOR="#EBEBEB" ROWSPAN=2> explanation </td>
        </tr>
        <tr>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> @p top </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> @p include_all </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> @p c.:fix_peptides </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 1 </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 2 </td>
            <td ALIGN="center" BGCOLOR="#EBEBEB"> sample 3 </td>
        </tr>
        <tr>
            <td ALIGN="center"> 0 </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> 1, 2, 3, 4 </td>
            <td ALIGN="center"> 2, 3, 4 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> all peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 1 </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> 1 </td>
            <td ALIGN="center"> 2 </td>
            <td ALIGN="center"> 1 </td>
            <td> single most abundant peptide </td>
        </tr>
        <tr>
            <td ALIGN="center"> 2 </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> 1, 2 </td>
            <td ALIGN="center"> 2, 3 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> two most abundant peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> 1, 2, 3 </td>
            <td ALIGN="center"> 2, 3, 4 </td>
            <td ALIGN="center"> (too few peptides) </td>
            <td> three most abundant peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> 1, 2, 3 </td>
            <td ALIGN="center"> 2, 3, 4 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> three or fewer most abundant peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 4 </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> 1, 2, 3, 4 </td>
            <td ALIGN="center"> (too few peptides) </td>
            <td ALIGN="center"> (too few peptides) </td>
            <td> four most abundant peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 4 </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> 1, 2, 3, 4 </td>
            <td ALIGN="center"> 2, 3, 4 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> four or fewer most abundant peptides </td>
        </tr>
        <tr>
            <td ALIGN="center"> 0 </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> 3 </td>
            <td> all peptides present in every sample </td>
        </tr>
        <tr>
            <td ALIGN="center"> 1 </td>
            <td ALIGN="center"> * </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> 3 </td>
            <td> single peptide present in most samples </td>
        </tr>
        <tr>
            <td ALIGN="center"> 2 </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 1, 3 </td>
            <td ALIGN="center"> (peptide 1 missing) </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> two peptides present in most samples </td>
        </tr>
        <tr>
            <td ALIGN="center"> 2 </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 1, 3 </td>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> two or fewer peptides present in most samples </td>
        </tr>
        <tr>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> no </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 1, 2, 3 </td>
            <td ALIGN="center"> (peptide 1 missing) </td>
            <td ALIGN="center"> (peptide 2 missing) </td>
            <td> three peptides present in most samples </td>
        </tr>
        <tr>
            <td ALIGN="center"> 3 </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> yes </td>
            <td ALIGN="center"> 1, 2, 3 </td>
            <td ALIGN="center"> 2, 3 </td>
            <td ALIGN="center"> 1, 3 </td>
            <td> three or fewer peptides present in most samples </td>
        </tr>
    </table>
</CENTER>

<B>Further considerations for parameter selection</B>

With @p best_charge_and_fractions and @p aggregate, there is a trade-off between comparability of protein abundances within a sample and of abundances for the same protein across different samples.\n
Setting @p best_charge_and_fraction may increase reproducibility between samples, but will distort the proportions of protein abundances within a sample. The reason is that ionization properties vary between peptides, but should remain constant across samples. Filtering by charge state can help to reduce the impact of feature detection differences between samples.\n
For @p aggregate, there is a qualitative difference between @p (intensity weighted) mean/median and @p sum in the effect that missing peptide abundances have (only if @p include_all is set or @p top is 0): @p (intensity weighted) mean and @p median ignore missing cases, averaging only present values. If low-abundant peptides are not detected in some samples, the computed protein abundances for those samples may thus be too optimistic. @p sum implicitly treats missing values as zero, so this problem does not occur and comparability across samples is ensured. However, with @p sum the total number of peptides ("summands") available for a protein may affect the abundances computed for it (depending on @p top), so results within a sample may become unproportional.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPProteinQuantifier :
  public TOPPBase
{
public:

  TOPPProteinQuantifier() :
    TOPPBase("ProteinQuantifier", "Compute peptide and protein abundances"),
    algo_params_(), proteins_(), peptides_(), columns_headers_(),
    spectral_counting_(false) {}

protected:

  typedef PeptideAndProteinQuant::PeptideQuant PeptideQuant;
  typedef PeptideAndProteinQuant::ProteinQuant ProteinQuant;
  typedef PeptideAndProteinQuant::SampleAbundances SampleAbundances;
  typedef PeptideAndProteinQuant::Statistics Statistics;
  typedef ProteinIdentification::ProteinGroup ProteinGroup;

  Param algo_params_; // parameters for PeptideAndProteinQuant algorithm
  ProteinIdentification proteins_; // protein inference results (proteins)
  vector<PeptideIdentification> peptides_; // protein inference res. (peptides)
  ConsensusMap::ColumnHeaders columns_headers_; // information about experimental design
  bool spectral_counting_; // quantification based on spectral counting?

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML,idXML"));
    registerInputFile_("protein_groups", "<file>", "", "Protein inference results for the identification runs that were used to annotate the input (e.g. via the ProteinInference tool).\nInformation about indistinguishable proteins will be used for protein quantification.", false);
    setValidFormats_("protein_groups", ListUtils::create<String>("idXML"));

    registerInputFile_("design", "<file>", "", "input file containing the experimental design", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    // output
    registerOutputFile_("out", "<file>", "", "Output file for protein abundances", false);
    setValidFormats_("out", ListUtils::create<String>("csv"));

    registerOutputFile_("peptide_out", "<file>", "", "Output file for peptide abundances", false);
    setValidFormats_("peptide_out", ListUtils::create<String>("csv"));

    registerOutputFile_("mztab", "<file>", "", "Output file (mzTab)", false);
    setValidFormats_("mztab", ListUtils::create<String>("mzTab"));

    // algorithm parameters:
    addEmptyLine_();
    Param temp = PeptideAndProteinQuant().getParameters();
    registerFullParam_(temp);

    registerStringOption_("greedy_group_resolution", "<choice>", "false", "Pre-process identifications with greedy resolution of shared peptides based on the protein group probabilities. (Only works with an idXML file given as protein_groups parameter).", false);
    setValidStrings_("greedy_group_resolution", ListUtils::create<String>("true,false"));
    registerFlag_("ratios", "Add the log2 ratios of the abundance values to the output. Format: log_2(x_0/x_0) <sep> log_2(x_1/x_0) <sep> log_2(x_2/x_0) ...", false);
    registerFlag_("ratiosSILAC", "Add the log2 ratios for a triple SILAC experiment to the output. Only applicable to consensus maps of exactly three sub-maps. Format: log_2(heavy/light) <sep> log_2(heavy/middle) <sep> log_2(middle/light)", false);
    registerTOPPSubsection_("format", "Output formatting options");
    registerStringOption_("format:separator", "<sep>", "", "Character(s) used to separate fields; by default, the 'tab' character is used", false);
    registerStringOption_("format:quoting", "<method>", "double", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
    setValidStrings_("format:quoting", ListUtils::create<String>("none,double,escape"));
    registerStringOption_("format:replacement", "<x>", "_", "If 'quoting' is 'none', used to replace occurrences of the separator in strings before writing", false);

    // registerSubsection_("algorithm", "Algorithm parameters section");
  }

  // Param getSubsectionDefaults_(const String& /* section */) const
  // {
  //     return PeptideAndProteinQuant().getParameters();
  // }


  /// Write a table of peptide results.
  void writePeptideTable_(SVOutStream& out, const PeptideQuant& quant, const ExperimentalDesign& ed)
  {
    // write header:
    out << "peptide" << "protein" << "n_proteins" << "charge";
    if (ed.getNumberOfSamples() <= 1)
    {
      out << "abundance";
    }
    else
    {
      for (Size i = 0; i < ed.getNumberOfSamples(); ++i)
      {
        out << "abundance_" + String(i+1);
      }
    }
    out << "fraction" << endl;

    bool best_charge_and_fraction = algo_params_.getValue("best_charge_and_fraction") == "true";
    for (auto const & q : quant) // loop over sequence->peptide data
    {
      if (q.second.total_abundances.empty())
      { 
        continue; // not quantified
      }
      StringList accessions;
      for (String acc : q.second.accessions)
      {
        accessions.push_back(acc.substitute('/', '_'));
      }
      String protein = ListUtils::concatenate(accessions, "/");

      if (best_charge_and_fraction)
      {
        // write individual abundances (one line for each charge state and fraction):
        for (auto const & fa : q.second.abundances)
        {
          const Size fraction = fa.first;
          for (auto const & ab : fa.second)
          {
            out << q.first.toString() << protein << accessions.size() << ab.first;

            for (size_t sample_id = 0; sample_id < ed.getNumberOfSamples(); ++sample_id)
            {
              // write abundance for the sample if it exists, 0 otherwise:
              SampleAbundances::const_iterator pos = ab.second.find(sample_id);
              out << (pos != ab.second.end() ? pos->second : 0.0);
            }
            out << fraction << endl; // output fraction
          }
        }
      }
      else
      {
        // write total abundances (accumulated over all charge states and fractions):
        out << q.first.toString() << protein << accessions.size() << 0;

        for (size_t sample_id = 0; sample_id < ed.getNumberOfSamples(); ++sample_id)
        {
          // write abundance for the sample if it exists, 0 otherwise:
          SampleAbundances::const_iterator pos = q.second.total_abundances.find(sample_id);
          out << (pos != q.second.total_abundances.end() ? pos->second : 0.0);
        }

        out << "all" << endl;
      }
    }
  }

  /// Write a table of protein results.
  void writeProteinTable_(SVOutStream& out, const ProteinQuant& quant, const ExperimentalDesign& ed)
  {
    const bool print_ratios = getFlag_("ratios");
    const bool print_SILACratios = getFlag_("ratiosSILAC");
    // write header:
    out << "protein" << "n_proteins" << "protein_score" << "n_peptides";
    if (ed.getNumberOfSamples() <= 1)
    {
      out << "abundance";
    }
    else
    {
      for (Size i = 0; i < ed.getNumberOfSamples(); ++i)
      {
        out << "abundance_" + String(i+1);
      }

      // TODO MULTIPLEXING: check if correct
      // if ratios-flag is set, print log2-ratios. ratio_1 <sep> ratio_x ....
      if (print_ratios)
      {
        for (Size i = 0; i < ed.getNumberOfSamples(); ++i)
        {
          out << "ratio_" + String(i+1);
        }
      }
      // if ratiosSILAC-flag is set, print SILAC log2-ratios, only if three
      if (print_SILACratios && ed.getNumberOfSamples() == 3)
      {
        for (Size i = 0; i < ed.getNumberOfSamples(); ++i)
        {
          out << "SILACratio_" + String(i+1);
        }
      }
    }

    out << endl;

    // mapping: accession of leader -> (accessions of grouped proteins, score)
    map<String, pair<StringList, double> > leader_to_group;
    if (!proteins_.getIndistinguishableProteins().empty())
    {
      for (auto group : proteins_.getIndistinguishableProteins()) //OMS_CODING_TEST_EXCLUDE
      {
        StringList& accessions = leader_to_group[group.accessions[0]].first;
        accessions = group.accessions;
        for (auto & acc : accessions)
        {
          acc.substitute('/', '_'); // to allow concatenation later
        }
        leader_to_group[group.accessions[0]].second = group.probability;
      }
    }

    for (auto const & q : quant)
    {
      if (q.second.total_abundances.empty())
      {
        continue; // not quantified
      }
      if (leader_to_group.empty())
      {
        out << q.first << 1;
        if (proteins_.getHits().empty())
        {
          out << 0;
        }
        else
        {
          vector<ProteinHit>::iterator pos = proteins_.findHit(q.first);
          out << pos->getScore();
        }
      }
      else
      {
        pair<StringList, double>& group = leader_to_group[q.first];
        out << ListUtils::concatenate(group.first, '/') << group.first.size()
            << group.second;
      }
      Size n_peptide = q.second.abundances.size();
      out << n_peptide;

      for (size_t sample_id = 0; sample_id < ed.getNumberOfSamples(); ++sample_id)
      {
        // write abundance for the sample if it exists, 0 otherwise:
        SampleAbundances::const_iterator pos = q.second.total_abundances.find(sample_id);
        out << (pos != q.second.total_abundances.end() ? pos->second : 0.0);
      }

      // if ratios-flag is set, print log2-ratios. ab1/ab0, ab2/ab0, ... , ab'n/ab0
      if (print_ratios)
      {
        double log2 = log(2.0);
        double ref_abundance = q.second.total_abundances.find(0)->second;
        out << 0; // =log(1)/log2;
        for (size_t sample_id = 1; sample_id < ed.getNumberOfSamples(); ++sample_id)

        {
          SampleAbundances::const_iterator pos = q.second.total_abundances.find(sample_id);
          out << (pos != q.second.total_abundances.end() ? log(pos->second / ref_abundance) / log2 : 0.0);
        }
      }
      // if ratiosSILAC-flag is set, print log2-SILACratios. Only if three maps are provided (triple SILAC).
      if (print_SILACratios && ed.getNumberOfSamples() == 3)
      {
        double light = q.second.total_abundances.find(0)->second;
        double middle = q.second.total_abundances.find(1)->second;
        double heavy = q.second.total_abundances.find(2)->second;

        double log2 = log(2.0);

        out << log(heavy / light) / log2
            << log(heavy / middle) / log2
            << log(middle / light) / log2;
      }
      out << endl;
    }
  }


  /// Write comment lines before a peptide/protein table.
  void writeComments_(SVOutStream& out, const ExperimentalDesign& ed, const bool proteins = true)
  {
    String what = (proteins ? "Protein" : "Peptide");
    bool old = out.modifyStrings(false);
    bool is_ibaq = algo_params_.getValue("method") == "iBAQ";
    out << "# " + what + " abundances computed from file '" +
      getStringOption_("in") + "'" << endl;
    StringList relevant_params;
    if (proteins) // parameters relevant only for protein output
    {
      relevant_params.push_back("method");
      if (!is_ibaq)
      {
        relevant_params.push_back("top:N");
        Size top = algo_params_.getValue("top:N");
        if (top != 1)
        {
          relevant_params.push_back("top:aggregate");
          if (top != 0)
          {
            relevant_params.push_back("top:include_all");
          }
        }
      }
    }
    relevant_params.push_back("best_charge_and_fraction"); // also for peptide output

    if (ed.getNumberOfSamples() > 1) // flags only for consensusXML input
    {
      relevant_params.push_back("consensus:normalize");
      if (proteins)
      {
        relevant_params.push_back("consensus:fix_peptides");
      }
    }

    String params;
    for (const String& str : relevant_params)
    {
      String value = algo_params_.getValue(str).toString();
      if (value != "false")
      {
        params += str + "=" + value + ", ";
      }
    }
    if (params.empty())
    {
      params = "(none)";
    }
    else
    {
      params.resize(params.size() - 2); // remove trailing ", "
    }
    out << "# Parameters (relevant only): " + params << endl;

    if (ed.getNumberOfSamples() > 1 && ed.getNumberOfLabels() == 1)
    {
      String desc = "# Files/samples associated with abundance values below: ";

      const auto& ms_section = ed.getMSFileSection();

      map<String, String> sample_id_to_filename;
      for (const auto& e : ms_section)
      {
        String ed_filename = File::basename(e.path);
        String ed_label = e.label;
        String ed_sample = e.sample;
        sample_id_to_filename[e.sample] = ed_filename; // should be 0,...,n_samples-1
      }

      for (Size i = 0; i < ed.getNumberOfSamples(); ++i)
      {
        if (i > 0)
        {
          desc += ", ";
        }
        desc += String(i + 1) + ": '" + sample_id_to_filename[String(i)] + "'";
      }
      out << desc << endl;
    }
    out.modifyStrings(old);
  }

  /// Write processing statistics.
  void writeStatistics_(const Statistics& stats)
  {
    OPENMS_LOG_INFO << "\nProcessing summary - number of...";
    if (spectral_counting_)
    {
      OPENMS_LOG_INFO << "\n...spectra: " << stats.total_features << " identified"
               << "\n...peptides: " << stats.quant_peptides
               << " identified and quantified (considering best hits only)";
    }
    else
    {
      OPENMS_LOG_INFO << "\n...features: " << stats.quant_features
               << " used for quantification, " << stats.total_features
               << " total (" << stats.blank_features << " no annotation, "
               << stats.ambig_features << " ambiguous annotation)"
               << "\n...peptides: "  << stats.quant_peptides
               << " quantified, " << stats.total_peptides
               << " identified (considering best hits only)";
    }
    if (!getStringOption_("out").empty())
    {
      bool include_all = algo_params_.getValue("top:include_all") == "true";
      Size top_n = algo_params_.getValue("top:N");
      OPENMS_LOG_INFO << "\n...proteins/protein groups: " << stats.quant_proteins
               << " quantified";
      if (top_n > 1)
      {
        if (include_all)
        {
          OPENMS_LOG_INFO << " (incl. ";
        }
        else
        {
          OPENMS_LOG_INFO << ", ";
        }
        OPENMS_LOG_INFO << stats.too_few_peptides << " with fewer than " << top_n
                 << " peptides";
        if (stats.n_samples > 1)
        {
          OPENMS_LOG_INFO << " in every sample";
        }
        if (include_all)
        {
          OPENMS_LOG_INFO << ")";
        }
      }
    }
    OPENMS_LOG_INFO << endl;
  }

  ExperimentalDesign getExperimentalDesignIds_(const String & design_file, const vector<ProteinIdentification> & proteins)
  {
    if (!design_file.empty()) // load experimental design file
    {
      return ExperimentalDesignFile::load(design_file, false);
      // TODO FRACTIONS: check if ed sane
    }
    else  // no design file provided
    {
      return ExperimentalDesign::fromIdentifications(proteins);
    }
  }

  ExperimentalDesign getExperimentalDesignFeatureMap_(const String & design_file, const FeatureMap & fm)
  {
    if (!design_file.empty()) // experimental design file
    {
      return ExperimentalDesignFile::load(design_file, false);
      // TODO FRACTIONS: check if ed sane
    }
    else  // no design given
    {
      return ExperimentalDesign::fromFeatureMap(fm);
    }
  }

  ExperimentalDesign getExperimentalDesignConsensusMap_(const String & design_file, const ConsensusMap & cm)
  {
    if (!design_file.empty()) // load experimental design file
    {
      return ExperimentalDesignFile::load(design_file, false);
      // TODO FRACTIONS: check if ed sane
    }
    else  // no design file provided
    {
      OPENMS_LOG_INFO << "No design file given. Trying to infer from consensus map." << std::endl;
      return ExperimentalDesign::fromConsensusMap(cm);
    }
  }

  ExitCodes main_(int, const char**) override
  {
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String peptide_out = getStringOption_("peptide_out");
    String mztab = getStringOption_("mztab");
    String design_file = getStringOption_("design");
    bool greedy_group_resolution = getStringOption_("greedy_group_resolution") == "true";

    if (out.empty() && peptide_out.empty())
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
                                                 OPENMS_PRETTY_FUNCTION,
                                                 "out/peptide_out");
    }

    String protein_groups = getStringOption_("protein_groups");
    if (!protein_groups.empty()) // read protein inference data
    {
      vector<ProteinIdentification> proteins;
      FileHandler().loadIdentifications(protein_groups, proteins, peptides_, {FileTypes::IDXML});
      if (proteins.empty() || 
          proteins[0].getIndistinguishableProteins().empty())
      {
        throw Exception::MissingInformation(
         __FILE__,
         __LINE__,
         OPENMS_PRETTY_FUNCTION,
         "No information on indistinguishable protein groups found in file '" + protein_groups + "'");
      }
      proteins_ = proteins[0]; // inference data is attached to first ID run
      if (greedy_group_resolution)
      {
        PeptideProteinResolution ppr{};
        ppr.buildGraph(proteins_, peptides_);
        ppr.resolveGraph(proteins_, peptides_);
      }
    }

    FileTypes::Type in_type = FileHandler::getType(in);

    PeptideAndProteinQuant quantifier;
    algo_params_ = quantifier.getParameters();
    Logger::LogStream nirvana; // avoid parameter update messages
    algo_params_.update(getParam_(), false, nirvana);
    // algo_params_.update(getParam_());
    quantifier.setParameters(algo_params_);

    // iBAQ works only with feature intensity values in consensusXML or featureXML files
    if (algo_params_.getValue("method") == "iBAQ" && in.hasSuffix("idXML"))
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "Invalid input: idXML for iBAQ, only consensusXML or featureXML are valid");
    }

    // iBAQ can only quantify proteins
    if (algo_params_.getValue("method") == "iBAQ" && !peptide_out.empty())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__,
                                        OPENMS_PRETTY_FUNCTION,
                                        "Invalid output: peptide_out can not be set when using iBAQ");
    }

    ExperimentalDesign ed;

    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap features;
      FileHandler().loadFeatures(in, features, {FileTypes::FEATUREXML});
      columns_headers_[0].filename = in;

      ed = getExperimentalDesignFeatureMap_(design_file, features);

      // protein inference results in the featureXML?
      if (protein_groups.empty() &&
          (features.getProteinIdentifications().size() == 1) &&
          (!features.getProteinIdentifications()[0].getHits().empty()))
      {
        proteins_ = features.getProteinIdentifications()[0];
      }
      quantifier.readQuantData(features, ed);
      quantifier.quantifyPeptides(peptides_); // quantify on peptide level
      quantifier.quantifyProteins(proteins_);
    }
    else if (in_type == FileTypes::IDXML)
    {
      spectral_counting_ = true;
      vector<ProteinIdentification> proteins;
      vector<PeptideIdentification> peptides;
      FileHandler().loadIdentifications(in, proteins, peptides, {FileTypes::IDXML});
      for (Size i = 0; i < proteins.size(); ++i)
      {
        columns_headers_[i].filename = proteins[i].getSearchEngine() + "_" + proteins[i].getDateTime().toString();
      }

      ed = getExperimentalDesignIds_(design_file, proteins);

      // protein inference results in the idXML?
      if (protein_groups.empty() && (proteins.size() == 1) && 
          (!proteins[0].getHits().empty()))
      {
        proteins_ = proteins[0];
      }
      quantifier.readQuantData(proteins, peptides, ed);
      quantifier.quantifyPeptides(peptides_); // quantify on peptide level
      quantifier.quantifyProteins(proteins_);
    }
    else // consensusXML
    {
      ConsensusMap consensus;
      FileHandler().loadConsensusFeatures(in, consensus, {FileTypes::CONSENSUSXML});
      columns_headers_ = consensus.getColumnHeaders();

      ed = getExperimentalDesignConsensusMap_(design_file, consensus);

      bool inference_in_cxml = false;
      // protein inference results in the consensusXML or from external ID-only file?
      if (protein_groups.empty() &&
          (consensus.getProteinIdentifications().size() == 1) &&
          consensus.getProteinIdentifications()[0].hasInferenceData())
      {
        proteins_ = consensus.getProteinIdentifications()[0];
        inference_in_cxml = true;
      }

      quantifier.readQuantData(consensus, ed);
      quantifier.quantifyPeptides(peptides_); // quantify on peptide level
      quantifier.quantifyProteins(proteins_);

      // write mzTab file
      if (!mztab.empty())
      {
        // annotate quants to protein(groups) for easier export in mzTab
        auto const & protein_quants = quantifier.getProteinResults();
        quantifier.annotateQuantificationsToProteins(protein_quants, proteins_);
        if (!inference_in_cxml)
        {
          auto& prots = consensus.getProteinIdentifications();
          prots.insert(prots.begin(), proteins_); // insert inference information as first protein identification
        }
        else
        {
          std::swap(consensus.getProteinIdentifications()[0], proteins_);
        }
        /*
        * TODO: maybe an assertion that the numbers of quantified proteins / ind. proteins match
        auto n_ind_prot = consensus.getProteinIdentifications()[0].getIndistinguishableProteins().size();
        cout << "MzTab Export: " << n_ind_prot << endl;
        */

        // fill MzTab with meta data and quants annotated in identification data structure
        const bool report_unmapped(true);
        const bool report_unidentified_features(false);
        const bool report_subfeatures(false);
        MzTabFile().store(mztab,
          consensus, 
          !inference_in_cxml,
          report_unidentified_features,
          report_unmapped,
          report_subfeatures);
      }
    }

    // output:
    String separator = getStringOption_("format:separator");
    String replacement = getStringOption_("format:replacement");
    String quoting = getStringOption_("format:quoting");
    if (separator.empty())
    {
      separator = "\t";
    }
    String::QuotingMethod quoting_method;
    if (quoting == "none")
    {
      quoting_method = String::NONE;
    }
    else if (quoting == "double")
    {
      quoting_method = String::DOUBLE;
    }
    else
    {
      quoting_method = String::ESCAPE;
    }
    if (!peptide_out.empty())
    {
      ofstream outstr(peptide_out.c_str());
      SVOutStream output(outstr, separator, replacement, quoting_method);
      writeComments_(output, ed, false);
      writePeptideTable_(output, quantifier.getPeptideResults(), ed);
      outstr.close();
    }
    if (!out.empty())
    {
      ofstream outstr(out.c_str());
      SVOutStream output(outstr, separator, replacement, quoting_method);
      writeComments_(output, ed);
      writeProteinTable_(output, quantifier.getProteinResults(), ed);
      outstr.close();
    }

    writeStatistics_(quantifier.getStatistics());

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPProteinQuantifier t;
  return t.main(argc, argv);
}

/// @endcond
