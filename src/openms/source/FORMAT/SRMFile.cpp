// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SRMFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <cstdio>

namespace OpenMS
{

std::vector<::OpenSwath::SwathMap> SRMFile::loadMzML(const String& file,
                                                   const String& /*tmp*/, 
                                                   std::shared_ptr<ExperimentalSettings>& exp_meta)
{
  OPENMS_LOG_INFO << "Loading SRM mzML " << file << std::endl;
  std::shared_ptr<PeakMap> full_exp(new PeakMap);
  FileHandler fh;
  // load full chromatograms
  fh.loadExperiment(file, *full_exp, {FileTypes::MZML});

  // Normalize / canonicalize chromatogram precursor/product m/z and nativeIDs
  for (Size i = 0; i < full_exp->getChromatograms().size(); ++i)
  {
    MSChromatogram & chrom = full_exp->getChromatograms()[i];
    double prec_mz = chrom.getPrecursor().getMZ();
    double prod_mz = chrom.getProduct().getMZ();

    // If precursor/product m/z are not set, try to read them from CV terms
    // (e.g. isolation window target m/z CV term MS:1000827) which some
    // mzML writers put into the <precursor>/<product> cvParams rather than
    // the Precursor/Product mz fields.
    if (!(prec_mz > 0))
    {
      const auto & prec_terms = chrom.getPrecursor().getCVTerms();
      auto it = prec_terms.find("MS:1000827");
      if (it != prec_terms.end() && !it->second.empty())
      {
        try
        {
          prec_mz = it->second[0].getValue().toString().toDouble();
          Precursor p = chrom.getPrecursor();
          p.setMZ(prec_mz);
          chrom.setPrecursor(p);
        }
        catch (...) {}
      }
    }

    if (!(prod_mz > 0))
    {
      const auto & prod_terms = chrom.getProduct().getCVTerms();
      auto it = prod_terms.find("MS:1000827");
      if (it != prod_terms.end() && !it->second.empty())
      {
        try
        {
          prod_mz = it->second[0].getValue().toString().toDouble();
          Product q = chrom.getProduct();
          q.setMZ(prod_mz);
          chrom.setProduct(q);
        }
        catch (...) {}
      }
    }

    // If we still don't have precursor/product m/z, fall back to parsing
    // them from the nativeID (handles legacy nativeID formats)
    if (!((prec_mz > 0) || (prod_mz > 0)))
    {
      String nid = chrom.getNativeID();
      double c_prec = -1.0, c_prod = -1.0;
      if (sscanf(nid.c_str(), "%*[^0-9+\\-]%lf%*[^0-9+\\-]%lf", &c_prec, &c_prod) >= 2)
      {
        Precursor p; p.setMZ(c_prec);
        Product q; q.setMZ(c_prod);
        chrom.setPrecursor(p);
        chrom.setProduct(q);
        prec_mz = c_prec;
        prod_mz = c_prod;
      }
    }

    // Ensure meta-values and a canonical nativeID exist for whichever values we have
    if (prec_mz > 0)
    {
      chrom.setMetaValue("precursor_mz", prec_mz);
    }
    if (prod_mz > 0)
    {
      chrom.setMetaValue("product_mz", prod_mz);
    }

    // Build a compact canonical nativeID depending on available values
    String nid_out;
    if (prec_mz > 0 && prod_mz > 0)
    {
      nid_out = String("precursor=") + String(prec_mz) + ",product=" + String(prod_mz);
    }
    else if (prod_mz > 0)
    {
      nid_out = String("product=") + String(prod_mz);
    }
    else if (prec_mz > 0)
    {
      nid_out = String("precursor=") + String(prec_mz);
    }

    if (!nid_out.empty()) chrom.setNativeID(nid_out);
  }

  exp_meta = full_exp;
  std::vector<::OpenSwath::SwathMap> swath_maps;
  ::OpenSwath::SwathMap sm;
  sm.sptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(full_exp);
  sm.lower = -1;
  sm.upper = -1;
  sm.center = -1;
  sm.imLower = -1;
  sm.imUpper = -1;
  sm.ms1 = false;
  swath_maps.push_back(sm);
  return swath_maps;
}

} // namespace OpenMS
