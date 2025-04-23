// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/RangeUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <map>

namespace OpenMS
{
  /**
    @brief Conversion class to convert chromatograms

    There are basically two methods implemented, conversion of
    chromatograms into spectra representation and vice versa.

    @ingroup Kernel
  */
  class ChromatogramTools
  {
public:
    /// @name Constructors and destructors
    //@{
    /// default constructor
    ChromatogramTools()
    {}

    /// copy constructor
    ChromatogramTools(const ChromatogramTools & /*rhs*/)
    {}

    /// destructor
    virtual ~ChromatogramTools()
    {}

    //@}

    /// @name Accessors
    //@{
    /**
      @brief converts the chromatogram to a list of spectra with instrument settings

      This conversion may be necessary as most of the spectra formats do not support
      chromatograms, except of mzML. However, most formats support e.g. SRM chromatogram
      as a list of spectra with instrument settings SRM and a separate spectrum for
      each data point. The disadvantage of storing chromatograms in spectra is its
      exhaustive memory consumption.
    */
    template <typename ExperimentType>
    void convertChromatogramsToSpectra(ExperimentType & exp)
    {
      for (std::vector<MSChromatogram >::const_iterator it = exp.getChromatograms().begin(); it != exp.getChromatograms().end(); ++it)
      {
        // for each peak add a new spectrum
        for (typename ExperimentType::ChromatogramType::const_iterator pit = it->begin(); pit != it->end(); ++pit)
        {
          typename ExperimentType::SpectrumType spec;

          // add precursor and product peaks to spectrum settings
          spec.getPrecursors().push_back(it->getPrecursor());
          spec.getProducts().push_back(it->getProduct());
          spec.setRT(pit->getRT());
          spec.setMSLevel(2);
          spec.setInstrumentSettings(it->getInstrumentSettings());
          spec.setAcquisitionInfo(it->getAcquisitionInfo());
          spec.setSourceFile(it->getSourceFile());

          // TODO implement others
          if (it->getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
          {
            spec.getInstrumentSettings().setScanMode(InstrumentSettings::SRM);
          }
          if (it->getChromatogramType() == ChromatogramSettings::SELECTED_ION_MONITORING_CHROMATOGRAM)
          {
            spec.getInstrumentSettings().setScanMode(InstrumentSettings::SIM);
          }

          // new spec contains one peak, with product m/z and intensity
          spec.emplace_back(it->getMZ(), pit->getIntensity());
          exp.addSpectrum(spec);
        }
      }

      exp.setChromatograms(std::vector<MSChromatogram >());
    }

    /**
      @brief converts e.g. SRM spectra to chromatograms

      This conversion is necessary to convert chromatograms, e.g. from SRM or MRM
      experiments to real chromatograms. mzML 1.1.0 has support for chromatograms
      which can be stored much more efficiently than spectra based chromatograms.
      However, most other file formats do not support chromatograms.

      @param exp the experiment to be converted.
      @param remove_spectra if set to true, the chromatogram spectra are removed from the experiment.
      @param force_conversion Convert even if ScanMode is not SRM or if there are no precursors (e.g. GC-MS data)
    */
    template <typename ExperimentType>
    void convertSpectraToChromatograms(ExperimentType & exp, bool remove_spectra = false, bool force_conversion = false)
    {
      typedef typename ExperimentType::SpectrumType SpectrumType;
      std::map<double, std::map<double, std::vector<SpectrumType> > > chroms;
      std::map<double, MSChromatogram > chroms_xic;
      for (typename ExperimentType::ConstIterator it = exp.begin(); it != exp.end(); ++it)
      {
        // TODO other types
        if (it->getInstrumentSettings().getScanMode() == InstrumentSettings::SRM || force_conversion)
        {
          // exactly one precursor and one product ion
          if (it->getPrecursors().size() == 1 && it->size() == 1)
          {
            chroms[it->getPrecursors().begin()->getMZ()][it->begin()->getMZ()].push_back(*it);
          }
          // Exactly one precursor and more than one product ion.
          // This is how some converters (e.g. ReAdW 4.0.2) store SRM data,
          // collecting all information from one precursor in a single
          // pseudo-spectrum
          else if (it->getPrecursors().size() == 1 && it->size() > 0)
          {
            for (Size peak_idx = 0; peak_idx < it->size(); peak_idx++)
            {
              // copy spectrum and delete all data, but keep metadata, then add single peak
              SpectrumType dummy = *it;
              dummy.clear(false);
              dummy.push_back((*it)[peak_idx]);
              chroms[it->getPrecursors().begin()->getMZ()][(*it)[peak_idx].getMZ()].push_back(dummy);
            }
          }
          // We have no precursor, so this may be a MS1 chromatogram scan (as encountered in GC-MS)
          else if (force_conversion)
          {
            for (auto& p : *it)
            {
              double mz = p.getMZ();
              ChromatogramPeak chr_p;
              chr_p.setRT(it->getRT());
              chr_p.setIntensity(p.getIntensity());
              if (chroms_xic.find(mz) == chroms_xic.end())
              {
                // new chromatogram
                chroms_xic[mz].getPrecursor().setMZ(mz);
                // chroms_xic[mz].setProduct(prod); // probably no product
                chroms_xic[mz].setInstrumentSettings(it->getInstrumentSettings());
                chroms_xic[mz].getPrecursor().setMetaValue("description", String("XIC @ " + String(mz)));
                chroms_xic[mz].setAcquisitionInfo(it->getAcquisitionInfo());
                chroms_xic[mz].setSourceFile(it->getSourceFile());
              }
              chroms_xic[mz].push_back(chr_p);
            }
          }
          else
          {
            OPENMS_LOG_WARN << "ChromatogramTools: need exactly one precursor (given " << it->getPrecursors().size() <<
            ") and one or more product ions (" << it->size() << "), skipping conversion of this spectrum to chromatogram. If this is a MS1 chromatogram, please force conversion (e.g. with -convert_to_chromatograms)." << std::endl;
          }
        }
        else
        {
          // This does not makes sense to warn here, because it would also warn on simple mass spectra...
          // TODO think what to to here
          //OPENMS_LOG_WARN << "ChromatogramTools: cannot convert other chromatogram spectra types than 'Selected Reaction Monitoring', skipping conversion." << std::endl;
          //
        }
      }

      // Add the XIC chromatograms
      for (auto & chrom: chroms_xic) exp.addChromatogram(chrom.second);

      // Add the SRM chromatograms
      typename std::map<double, std::map<double, std::vector<SpectrumType> > >::const_iterator it1 = chroms.begin();
      for (; it1 != chroms.end(); ++it1)
      {
        typename std::map<double, std::vector<SpectrumType> >::const_iterator it2 = it1->second.begin();
        for (; it2 != it1->second.end(); ++it2)
        {
          typename ExperimentType::ChromatogramType chrom;
          chrom.setPrecursor(*it2->second.begin()->getPrecursors().begin());
          Product prod;
          prod.setMZ(it2->first);
          chrom.setProduct(prod);
          chrom.setInstrumentSettings(it2->second.begin()->getInstrumentSettings());
          chrom.setAcquisitionInfo(it2->second.begin()->getAcquisitionInfo());
          chrom.setSourceFile(it2->second.begin()->getSourceFile());

          typename std::vector<SpectrumType>::const_iterator it3 = it2->second.begin();
          for (; it3 != it2->second.end(); ++it3)
          {
            typename ExperimentType::ChromatogramType::PeakType p;
            p.setRT(it3->getRT());
            p.setIntensity(it3->begin()->getIntensity());
            chrom.push_back(p);
          }

          chrom.setNativeID("chromatogram=" + it2->second.begin()->getNativeID());               // TODO native id?
          chrom.setChromatogramType(ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM);
          exp.addChromatogram(chrom);
        }
      }

      if (remove_spectra)
      {
        exp.getSpectra().erase(remove_if(exp.begin(), exp.end(), HasScanMode<SpectrumType>(InstrumentSettings::SRM)), exp.end());
      }
    }

    //@}
  };
} // namespace OpenMS

