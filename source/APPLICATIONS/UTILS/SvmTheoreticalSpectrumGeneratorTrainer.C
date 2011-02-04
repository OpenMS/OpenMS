// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SvmTheoreticalSpectrumGeneratorTrainer.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>

//#define SGNT_DEBUG
#undef SGNT_DEBUG

using namespace std;
using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_SpectrumGeneratorSvmTheoreticalSpectrumGeneratorTrainer SpectrumGeneratorSvmTheoreticalSpectrumGeneratorTrainer

  @brief Trainer for Svm model as input for SvmSpectrumGenerator.

  This application requires a list of annotated spectra and generates a trained svm for each selected
  ion type (i.e. a,b,c,x,y,z + losses). It uses the libsvm library and all svmlib parameters are accessible
  by user. Please refer to the libsvm manuals for more detailed explanation of the parameters. The default
  values are choses as in the svm-training tool delivered with libsvm.

  @note This tool is experimental!

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SpectrumGeneratorNetworkTrainer.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class SvmTheoreticalSpectrumGeneratorTrainerTOPP
  : public TOPPBase
{
  typedef SvmTheoreticalSpectrumGenerator::IonType IonType;

public:
  SvmTheoreticalSpectrumGeneratorTrainerTOPP() :
      TOPPBase("SvmTheoreticalSpectrumGeneratorTrainerTOPP", "Trainer for SVM models as input for SvmSpectrumGenerator", false)
  {
  }

    protected:
  void registerOptionsAndFlags_()
  {
    // I/O settings
    registerInputFile_("in_spectra", "<file>", "", "Input Training Spectra in mzML", true);
    registerInputFile_("in_identifications", "<file>", "", "Input file with corresponding sequences in IdXML", true);
    registerOutputFile_("model_output_file", "<file>", "",
                         "Name for output files. For each ion_type one file <filename>_residue_loss_charge.svm and one <filename>.info which has to be passed to the SvmTheoretical SpectrumGenerator", true);
    registerIntOption_("precursor_charge", "<Int>", 2, "Precursor charge state used for model training", false);
    setMinInt_("precursor_charge",1);
    setMaxInt_("precursor_charge",3);

    registerSubsection_("algorithm", "");
  }

  Param getSubsectionDefaults_(const String& /* section*/) const
  {
    Param tmp = SvmTheoreticalSpectrumGeneratorTrainer().getDefaults();
    return tmp;
  }



  ExitCodes main_(int , const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in_spectra = getStringOption_("in_spectra");
    String in_identifications = getStringOption_("in_identifications");
    String outfile = getStringOption_("model_output_file");
    Size precursor_charge = getIntOption_("precursor_charge");

    //-------------------------------------------------------------
    // init SvmTheoreticalSpectrumGeneratorTrainer
    //-------------------------------------------------------------
    SvmTheoreticalSpectrumGeneratorTrainer trainer;

    Param param = getParam_().copy("algorithm:",true);
    trainer.setParameters(param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    PeakMap map;
    MzMLFile().load(in_spectra,map);

    std::vector<PeptideIdentification>pep_ids;
    std::vector<ProteinIdentification> prot_ids;
    String tmp_str;
    IdXMLFile().load(in_identifications,prot_ids, pep_ids, tmp_str);

    IDMapper idmapper;
    Param par;
    par.setValue("rt_tolerance", 0.001);
    par.setValue("mz_tolerance", 0.001);
    idmapper.setParameters(par);
    idmapper.annotate(map, pep_ids, prot_ids);

    //generate vector of annotations
    std::vector<AASequence>annotations;
    PeakMap::iterator it;
    for(it=map.begin(); it!=map.end(); ++it)
    {
      annotations.push_back(it->getPeptideIdentifications()[0].getHits()[0].getSequence());
    }

    trainer.trainModel(map, annotations, outfile, precursor_charge);
    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  SvmTheoreticalSpectrumGeneratorTrainerTOPP tool;
  return tool.main(argc, argv);
}

/// @endcond

