// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff, Alexander Haupt$
// --------------------------------------------------------------------------


#include <OpenMS/SIMULATION/LCMSSample.h>

using namespace std;


namespace OpenMS
{
  //-------------------------------------------------------------
  // Constructor/Deconstructor
  //-------------------------------------------------------------

  LCMSSample::LCMSSample()
  :	DefaultParamHandler("LCMSSample"),
    min_detect_(0.0)
  {
    defaults_.setValue("missed_cleavages",0,"maximum number of missed cleavages");
    defaults_.setValue("min_peptide_length",0,"minimum peptide length after digestion");
    defaults_.setValue("min_detect",0.5,"minimum peptide detectability accepted");

		defaultsToParam_();
  }

  LCMSSample::~LCMSSample() { }


  LCMSSample::LCMSSample(const LCMSSample& source)
    : DefaultParamHandler(source), peptides_(source.peptides_),
      proteins_(source.proteins_)
  {
    setParameters( source.getParameters() );
    updateMembers_();
  }

  LCMSSample& LCMSSample::operator = (const LCMSSample& source)
  {
    if (&source == this) return *this;

    setParameters( source.getParameters() );
    updateMembers_();

    peptides_ = source.peptides_;
    proteins_ = source.proteins_;

    return *this;
  }

  // Load proteins from FASTA file
  void LCMSSample::loadFASTA(const String filename)
  {
    FASTAFile fastafile;
    typedef std::vector< FASTAFile::FASTAEntry > FASTAdata;
    FASTAdata fastadata;

    // load FASTA file contents
    fastafile.load(filename, fastadata);
    // reserve vector space
    proteins_.reserve(fastadata.size());

    // add data from file to protein storage
    String::size_type index;
    int relativeQuantity;
    for (FASTAdata::iterator it = fastadata.begin(); it != fastadata.end(); ++it)
    {
      // remove all ambiguous characters from FASTA entry
      it->sequence.remove('X');
      it->sequence.remove('B');
      it->sequence.remove('Z');

      // Look for a relative quantity given in the first line of a FASTA entry
      index = (it->identifier).find_first_of("#");
      // if found, extract and set relative quantity accordingly
      if (index != string::npos)
      {
        stringstream strm ((it->identifier).substr(0, index));
        strm >> relativeQuantity;
      }
      else
      {
        relativeQuantity = 1;
      }

      proteins_.push_back(make_pair(it->sequence, relativeQuantity));
    }

    cout << endl;
    cout << "Done." << flush;
    cout << " (" << fastadata.size() << " proteins loaded)";

    cout << endl;

  }

  void LCMSSample::digest()
  {
    EnzymaticDigestion digestion;
    digestion.setEnzyme(EnzymaticDigestion::TRYPSIN);

    UInt missedCleavages  = param_.getValue("missed_cleavages");
    UInt minPeptideLength = param_.getValue("min_peptide_length");

    digestion.setMissedCleavages( missedCleavages );

    // initialize support vector machine
    LibSVMEncoder encoder;
    svm_problem* training_data = NULL;
    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;

		if (File::readable(DtModelFile_))
		{
    	//cout << "Loading svm model: " << DtModelFile_ << endl;
    	svm_.loadModel(DtModelFile_);
			//cout << "Done. " << endl;
    }
		
		// load additional parameters
		if (svm_.getIntParameter(KERNEL_TYPE) == OLIGO)
    {
      String add_paramfile = DtModelFile_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        cout << "SVM parameter file " << DtModelFile_ << " not found or not readable" << endl;
        cout << "Not performing detectability prediction!" << endl;
      }

      Param additional_parameters;
      additional_parameters.load(add_paramfile);

      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm_.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
         cout << "No border length defined in additional parameters file." << endl;
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm_.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
         cout << "No k-mer length defined in additional parameters file. Aborting detectability prediction!" << endl;
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();

      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm_.getIntParameter(KERNEL_TYPE) == OLIGO)
      {
        cout << "No sigma defined in additional parameters file. Aborting detectability prediction!" << endl;
      }

      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }

		if (File::readable(DtModelFile_))
		{
    	svm_.setParameter(BORDER_LENGTH, (Int) border_length);
    	svm_.setParameter(SIGMA, sigma);
    	// to obtain probabilities
    	svm_.setParameter(PROBABILITY, 1);
		}
      // loading training data
    String sample_file = DtModelFile_ + "_samples";
    if (File::readable(sample_file))
    {
    	training_data = encoder.loadLibSVMProblem(sample_file);
    	svm_.setTrainingSample(training_data);
    }
    cout << "Digesting...        " << endl;

    ProgressLogger plog;
    plog.setLogType(ProgressLogger::CMD);
    plog.startProgress(0, proteins_.size(), "Protein digestion");
    UInt c = 0;

    vector<String> all_peptides;
    vector<String> filtered_peptides;

    vector<AASequence> digestionProducts;
    stringstream strm;

    // Iterate through proteins and digest them
    for (SampleProteinsConstIt protein = proteins_.begin();
         protein != proteins_.end();
         ++protein)
    {
      AASequence sequence = AASequence(protein->first);
      digestion.digest(sequence, digestionProducts);

      for (vector<AASequence>::const_iterator dp_it = digestionProducts.begin();
           dp_it != digestionProducts.end();
           ++dp_it)
      {
        if (dp_it->size() < minPeptideLength) continue;

        strm << *dp_it;
        all_peptides.push_back( strm.str() );
        strm.str("");
      }

      // remove peptides with detectability likelihood below threshold
      filterForDetectability_(all_peptides, filtered_peptides, k_mer_length);
      //cout << "Kept " << filtered_peptides.size() << " out of " << all_peptides.size()  << endl;

      for (vector<String>::iterator it = filtered_peptides.begin();
            it != filtered_peptides.end();
            ++it)
      {
        peptides_[ *it ] += protein->second;
      }

      all_peptides.clear();
      digestionProducts.clear();
      filtered_peptides.clear();

      plog.setProgress(++c);
    }
  #ifdef DEBUG_SIM
    cout << "Done." << endl;
    cout << "Kept " << peptides_.size() << " peptides in total." << endl;
    cout << endl;
  #endif

    delete training_data;
    plog.endProgress();
  }

  void LCMSSample::filterForDetectability_(const vector<String>& all_peptides, vector<String>& filtered_peptides, UInt k_mer_length)
  {

    if (!File::readable(DtModelFile_))
    {
      cout << "Peptide detectability prediction disabled." << endl;
      cout << "All peptides will appear in LC-MS map. " << endl;
      filtered_peptides = all_peptides;
      return;
    }

    //cout << "Predicting peptide detectabilities..    " << endl;

    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";

    LibSVMEncoder encoder;

    // Encoding test data
    vector<DoubleReal> probs;
    probs.resize(all_peptides.size(), 0);

    svm_problem* prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(all_peptides, probs,
                                                                                     k_mer_length,
                                                                                     allowed_amino_acid_characters,
                                                                                     svm_.getIntParameter(BORDER_LENGTH));

    vector<DoubleReal> labels;
    vector<DoubleReal> detectabilities;
    svm_.getSVCProbabilities(prediction_data, detectabilities, labels);

    //cout << "Done." << endl;

    delete prediction_data;

    #ifdef DEBUG_SIM
    cout << "----------------------------------------------------------------" << endl;
    cout << "Predicted detectabilities:" << endl;
    #endif

    UInt c = 0;
    for (vector<String>::const_iterator it = all_peptides.begin();
           it != all_peptides.end();
           ++it)
    {
      #ifdef DEBUG_SIM
      cout << detectabilities[c] << " " << min_detect_ << endl;
      #endif
      if (detectabilities[c] > min_detect_) filtered_peptides.push_back(*it);
      ++c;
    }

    #ifdef DEBUG_SIM
    cout << "----------------------------------------------------------------" << endl;
    #endif
  }

  void LCMSSample::updateMembers_()
  {
    min_detect_ = param_.getValue("min_detect");
  }

  //-------------------------------------------------------------
  // Data output
  //-------------------------------------------------------------

  void LCMSSample::printProteins() const
  {
    UInt counter = 0;

    for(SampleProteinsConstIt it = proteins_.begin(); it != proteins_.end(); ++it)
    {
      cout << it->first << endl;
      counter++;
    }
    cout << counter << " protein(s) in the sample" << endl;
  }

  void LCMSSample::printPeptides() const
  {
    UInt counter = 0;

    for(ConstIterator it = peptides_.begin(); it != peptides_.end(); ++it)
    {
      cout << it->first << " (" << it->second << ")" << endl;
      counter += it->second;
    }
    cout << counter << " peptides in the sample" << endl;
  }

}
