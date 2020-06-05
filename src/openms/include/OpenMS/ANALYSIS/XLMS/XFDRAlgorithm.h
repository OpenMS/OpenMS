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
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>

namespace OpenMS
{

  //-------------------------------------------------------------
  // Doxygen docu
  //-------------------------------------------------------------

  /**
      @brief Calculates false discovery rate estimates on crosslink identifications.

      This tool calculates and FDR estimate for crosslink identifications, which are produced by OpenPepXL.
      The method employed currently is identical to the target-decoy approach used by xProphet (Walzthoeni et al., 2012).
      Consequently, this tool can also consume xquest.xml files (produced either by OpenPepXL or xQuest). The tool supports
      output in the idXML and mzIdentML formats.

      @experimental This tool is work in progress and usage and input requirements might change.


  */

  class OPENMS_DLLAPI XFDRAlgorithm :
  public DefaultParamHandler, public ProgressLogger
  {

  public:

    /// Exit codes
    enum ExitCodes
    {
      EXECUTION_OK,
      ILLEGAL_PARAMETERS,
      UNEXPECTED_RESULT
    };

    /// Default constructor
    XFDRAlgorithm();

    /// Default destructor
    ~XFDRAlgorithm() override;

    /**
     * @brief Performs the main function of this class, the FDR estimation for cross-linked peptide experiments

     * @param peptide_ids The PeptideIdentifications from an XL-MS experiment

     */
    ExitCodes run(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id);

    /**
    * @brief Checks whether the parameters of the object are valid
    *  @return ExitCode EXECUTION_OK if they are valid, ILLEGAL_PARAMETERS otherwise
    */
    ExitCodes validateClassArguments() const;

private:
    void updateMembers_() override;

    /**
   * @brief Prepares vector of PeptideIdentification such that it can be processed downstream.
   * The encompassed steps are:
   *  * Set min_score_ and max_score_ encountered in the data
   *  * Ensure that crosslink_type and crosslink_rank are available in the PeptideIdentification
   *  * Define the crosslink as either inter/or intraprotein
   *  * Set the identifier of the Peptide Identification if there is only one protein identification
   *
   */
    void initDataStructures_(std::vector<PeptideIdentification>& peptide_ids, ProteinIdentification& protein_id);

    /**
     * @brief Inspects PeptideIdentification pep_id and assigns all cross-link types that this identification belongs to
     * @param pep_id Peptide ID to be assigned.
     * @param types Result vector containing the names of the crosslink classes
     */
    static void assignTypes_(PeptideHit& ph, StringList& types);

    /** Target counting as performed by the xProphet software package
     *
     * @brief xprophet  method for target hits counting as implemented in xProphet
     * @param cum_histograms Cumulative score distributions
     */
    void fdr_xprophet_(std::map< String, Math::Histogram<> >& cum_histograms,
                      const String& targetclass, const String& decoyclass, const String& fulldecoyclass,
                      std::vector< double >& fdr, bool mono);

    /**
    * @brief Calculates the qFDR values for the provided FDR values, assuming that the FDRs are sorted by score in the input vector
    * @param fdr Vector with FDR values which should be used for qFDR calculation
    * @param qfdr Result qFDR values
    */
    static void calc_qfdr_(const std::vector< double >& fdr, std::vector< double >& qfdr);

    void findTopUniqueHits_(std::vector<PeptideIdentification>& peptide_ids);

    void writeArgumentsLog_() const;

    String getId_(const PeptideHit& ph) const;

    static Size getMinIonsMatched_(const PeptideHit& ph)
    {
      Size alpha_ions = Size(ph.getMetaValue("matched_linear_alpha")) + Size(ph.getMetaValue("matched_xlink_alpha"));
      Size beta_ions = Size(ph.getMetaValue("matched_linear_beta")) + Size(ph.getMetaValue("matched_xlink_beta"));
      return std::min(alpha_ions, beta_ions);
    }

    inline static void setIntraProtein_(PeptideHit& ph, const bool value)
    {
      ph.setMetaValue("XFDR:is_intraprotein", DataValue(value ? "true" : "false"));
    }

    inline static void setInterProtein_(PeptideHit& ph, const bool value)
    {
      ph.setMetaValue("XFDR:is_interprotein", DataValue(value ? "true" : "false"));
    }

    /**
     *  @brief Determines whether the Petide Evidences belong to the same protein, modulo decoy
     */
    static bool isSameProtein_(
            String prot1,
            String prot2,
            const String &decoy_string)
    {
      prot1.substitute(decoy_string, "");
      prot2.substitute(decoy_string, "");
      assert( ! prot1.hasSubstring(decoy_string));
      assert( ! prot2.hasSubstring(decoy_string));
      return prot1 == prot2;
    }

    // Score range for this of the tool
    Int min_score_;
    Int max_score_;

    // unique top hits
    std::vector<String> unique_ids_;
    std::vector<double> unique_id_scores_;

    // maps index of peptide id all_pep_ids_ to vector of cross link class
    std::map<String, std::vector<String>> cross_link_classes_;

    // Program arguments
    String decoy_string_;
    double arg_mindeltas_;
    double arg_minborder_;
    double arg_maxborder_;
    Int arg_minionsmatched_;
    double arg_minscore_;
    bool arg_uniquex_;
    bool arg_no_qvalues_;
    double arg_binsize_;

    // Names of the class parameters
    static const String param_decoy_string_;
    static const String param_minborder_;
    static const String param_maxborder_;
    static const String param_mindeltas_;
    static const String param_minionsmatched_;
    static const String param_uniquexl_;
    static const String param_no_qvalues_;
    static const String param_minscore_;
    static const String param_binsize_;

    // Constants related to particular crosslink classes
    static const String crosslink_class_intradecoys_;
    static const String crosslink_class_fulldecoysintralinks_;
    static const String crosslink_class_interdecoys_;
    static const String crosslink_class_fulldecoysinterlinks_;
    static const String crosslink_class_monodecoys_;
    static const String crosslink_class_intralinks_;
    static const String crosslink_class_interlinks_;
    static const String crosslink_class_monolinks_;
    static const String crosslink_class_decoys_;
    static const String crosslink_class_targets_;
    static const String crosslink_class_hybriddecoysintralinks_;
    static const String crosslink_class_hybriddecoysinterlinks_;
  };
}
