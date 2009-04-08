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
// $Authors: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

// General ideas:
//   * aggressively split the run method into several sub-units
//   * move as many parts as possible into separate classes
//     so we can reuse it for e.g. MALDI-MS Sim
//

#include <OpenMS/SIMULATION/LCMSSim.h>

using namespace std;

namespace OpenMS{

  LCMSSim::LCMSSim()
  : DefaultParamHandler("LCMSSim"),
  	ProgressLogger(),
    sample_(), RTModelFile_(), exp_(), gradientTime_(0.0), rt_sampling_(0.0),
    msAccuracy_(0), msBinSize_(0.0), mzMeanError_(0), mzStdDevError_(0), intMeanError_(0),
    intStdDevError_(0), maxMapMZ_(0.0), minMapMZ_(0.0), mean_scaling_(0.0), ion_count_(0),
    allowedMods_(), current_feature_(), features_(),distortion_(0.0),symmetry_up_(-100.0),symmetry_down_(100.0),
    changed_scans_(), allow_overlaps_(1)
  {

    // Column settings
    defaults_.setValue("total_gradient_time",2000.0,"the duration (in seconds) of the gradient");
    defaults_.setValue("rt_sampling",2.0,"Time interval [s] between consecutive scans");

    // MS instrument settings
    defaults_.setValue("ion_model",1,"Ionization model (0=SIMPLE, 1=ESI)");
    defaults_.setValue("peak_fwhm",0.5,"FWHM (full width at half maximum) of simulated peaks (Da).");
    defaults_.setValue("mz_sampling",0.12,"MS hardware resolution (e.g. bin size in m/z).");
  // 	defaults_.setValue("mz_accuracy",0.001,"MS hardware accuracy (ppm)");

    defaults_.setValue("mz_upperlimit",2500.0,"Upper m/z detecter limit.");
    defaults_.setValue("mz_lowerlimit",200.0,"Lower m/z detecter limit.");

    defaults_.setValue("allow_overlaps",1,"Disables simulation of overlapping signals.");

    // noise params
    // m/z error
    defaults_.setValue("mz_error_mean",0.0,"Average systematic m/z error (Da)");
    defaults_.setValue("mz_error_stddev",0.0,"Standard deviation for m/z errors. Set to 0 to disable simulation of m/z errors.");
    // intensity error
    defaults_.setValue("int_error_mean",0,"Average systematic intensity error.");
    defaults_.setValue("int_error_stddev",0.0,"Standard deviation for peak intensities (relative to peak height). Set to 0 to disable intensity errors.");
    // shot noise
    defaults_.setValue("noise_rate",0,"Poisson rate of shot noise. Set to 0 to disable simulation of shot noise.");
    defaults_.setValue("noise_int_mean",50.0,"Shot noise intensity mean.");
    // rt error
    defaults_.setValue("rt_shift_mean",0,"Mean shift in retention time [s]");
    defaults_.setValue("rt_shift_stddev",50,"Standard deviation of shift in retention time [s]");
    // contaminants
    defaults_.setValue("perc_contaminants",10.0,"Percentage of chemical contaminants");
    // posttranslational modifications
    defaults_.setValue("ptms_on",0,"Modeling of posttranslational modifications (0 = disabled) ");
    // baseline
    defaults_.setValue("baseline_scaling",0.0,"Scale of baseline (zero disables baseline simulation)");
    // column condition
    defaults_.setValue("column_condition","medium","LC condition: good | medium | poor");

    // random seed (advanced paramter?)
    defaults_.setValue("random_seed",0,"Number used to initialize the random number generator, for reproducible results (0 = time used as random seed)" /*,true*/);

    defaultsToParam_();

    // initialize random generator
    gsl_rng_default_seed = time(0);
    rand_gen_ = gsl_rng_alloc(gsl_rng_mt19937);

    // read chemical modifications
    readFromModFile_();

  }

  LCMSSim::~LCMSSim()
  {
    gsl_rng_free(rand_gen_);
  }

  LCMSSim::LCMSSim(const LCMSSim& source)
    : DefaultParamHandler(source),
  		ProgressLogger(source)
  {
		// initialize random generator
    gsl_rng_default_seed = time(0);
    rand_gen_ = gsl_rng_alloc(gsl_rng_mt19937);

    setParameters( source.getParameters() );
    updateMembers_();
  }

  LCMSSim& LCMSSim::operator = (const LCMSSim& source)
  {
    if (&source == this) return *this;

		// initialize random generator
    gsl_rng_default_seed = time(0);
    rand_gen_ = gsl_rng_alloc(gsl_rng_mt19937);
    setParameters( source.getParameters() );
    updateMembers_();

    return *this;
  }


  void LCMSSim::readFromModFile_()
  {
    std::ifstream in("ptms.dat");
    if (!in)
    {
      cout << "Could not read file <ptms.dat>." << endl;
      cout << "No modifications are modelled." << endl;
      return;
    }

    // Format: Name, Residue(s), Rel. Abundance, Formula, shift (space separated)
    vector<String> parts(5);
    String line;
    char delimiter = ' ';

    while (getline(in,line,'\n'))
    {
      line.trim();

      // skip comment line
      if ( line.empty() ||  line.hasPrefix("#") ) continue;

      if (line.has('\t'))
      {
        delimiter = '\t';
      }
      else
      {
        delimiter = ' ';
      }

      line.split(delimiter,parts);

      String res = parts[1];
      vector<String> residues;
      res.split(',',residues);

      PTM ptm;
      ptm.name_         = parts[0];
      ptm.abundance_ = parts[2].toDouble();
      ptm.formula_      = parts[3];

      // positive or negative shift?
      if (parts[4] == "+")
      {
        ptm.shift_ = true;
      }
      else
      {
        ptm.shift_ = false;
      }

      for (UInt i=0;i<residues.size();++i)
      {
        allowedMods_.insert(make_pair(residues[i],ptm));
      }
    }
  }

  void LCMSSim::readFromContaminationFile_(std::vector<EmpiricalFormula> & vef)
  {
    ifstream in("contaminants.dat");
    if (! File::readable("contaminants.dat") )
    {
      cout << "Could not read file <contaminants.dat>." << endl;
      cout << "No contaminations are modelled." << endl;
      return;
    }

    // Format: Name, Residue(s) , Rel. Abundance, Formula (tab separated)
    vector<String> strings(2);
    String line;
    char delimiter = '\t'; // Metabolite names can contain spaces, so we only allow tab-separation'

    while (getline(in,line,'\n'))
    {
      line.trim();

      if ( line.empty() ||  line.hasPrefix("#") )
      {
        continue;
      }

      line.split(delimiter,strings);

      if (strings.size() != 2)
      {
        cout << "Parse error in " << line << endl;
        continue;
      }

      try
      {
        vef.push_back(EmpiricalFormula(strings[1]));
      }
      catch (Exception::ParseError ex)
      {
        cout << "Parse error (" << ex.getMessage();
        cout << ") in formula: " << strings[1] << endl;
        cout << "Skipping this entry. " << endl;
      }
    }
  }

  void LCMSSim::predictRT_(RTTable& rtTable)
  {
    if (RTModelFile_ == "none")
    {
      // User didn't set a RT model file, so we assume he/she doesn't
      // want RT prediction at all. All peptides are going to elute
      // in the middile of the map.
      CoordinateType mid_rt = 0.5;

      cout << "RT prediction disabled." << endl;
      cout << "All peptides eluting at: " << mid_rt << endl;

      for (PeptideSequences::const_iterator seq_it = sample_.getPeptideSequences().begin();
            seq_it != sample_.getPeptideSequences().end();
            ++seq_it)
      {
        rtTable.insert( make_pair(mid_rt, seq_it) );
      }

      return;
    }
    else if(RTModelFile_ == "1D")
    {
      CoordinateType no_rt = 0.0;

      cout << "LC step was disabled." << endl;
      cout << "All peptides eluting at instantly " << endl;

      for (PeptideSequences::const_iterator seq_it = sample_.getPeptideSequences().begin();
            seq_it != sample_.getPeptideSequences().end();
            ++seq_it)
      {
        rtTable.insert( make_pair(no_rt, seq_it) );
      }

      return;
    }


    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    SVMWrapper svm;
    LibSVMEncoder encoder;
    vector<DoubleReal> predicted_retention_times;

    svm_problem* training_data = NULL;
    svm_problem* prediction_data = NULL;

    UInt k_mer_length = 0;
    DoubleReal sigma = 0.0;
    UInt border_length = 0;

    cout << "Predicting RT..    " << endl;

    // not that elegant...
    vector< String > peptidesVector;

    for(PeptideSequences::const_iterator seq_it = sample_.getPeptideSequences().begin();
         seq_it != sample_.getPeptideSequences().end();
         ++seq_it)
    {
      peptidesVector.push_back(seq_it->first);
    }

    svm.loadModel(RTModelFile_);

     // load additional parameters
    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
    {
      String add_paramfile = RTModelFile_ + "_additional_parameters";
      if (! File::readable( add_paramfile ) )
      {
        cout << "SVM parameter file " << add_paramfile << " not found or not readable" << endl;
        cout << "Aborting RT prediction!" << endl;
        //TODO where is the "return;" ?
        //TODO why is rtTable not filled with default values first?! (empty rtTable leads to empty map?! -->  see LCMSSim::run())
				// Ole: these are valid suggestions, but they should be adressed when the simulator is modified to 
				// incorporate MS/MS spectra.
      }

      Param additional_parameters;
      additional_parameters.load(add_paramfile);

      if (additional_parameters.getValue("border_length") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
         cout << "No border length defined in additional parameters file. Aborting RT prediction!" << endl;
         return;
      }
      border_length = ((String)additional_parameters.getValue("border_length")).toInt();
      if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        cout << "No k-mer length defined in additional parameters file. Aborting RT prediction!" << endl;
        return;
      }
      k_mer_length = ((String)additional_parameters.getValue("k_mer_length")).toInt();

      if (additional_parameters.getValue("sigma") == DataValue::EMPTY
          && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        cout << "No sigma defined in additional parameters file. Aborting RT prediction!" << endl;
        return;
      }

      sigma = ((String)additional_parameters.getValue("sigma")).toFloat();
    }

    svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) border_length);
    svm.setParameter(SVMWrapper::SIGMA, sigma);

    // Encoding test data
    vector<DoubleReal> rts;
    rts.resize(peptidesVector.size(), 0);
    prediction_data = encoder.encodeLibSVMProblemWithOligoBorderVectors(peptidesVector, rts, k_mer_length, allowed_amino_acid_characters, border_length);

    // loading training data
    String sample_file = RTModelFile_ + "_samples";
    if (! File::readable( sample_file ) )
    {
      cout << "SVM sample file " << sample_file << " not found or not readable" << endl;
      cout << "Aborting RT prediction!" << endl;
    }
    training_data = encoder.loadLibSVMProblem(sample_file);

    svm.setTrainingSample(training_data);
    svm.predict(prediction_data, predicted_retention_times);

    cout << "Done." << endl;
    delete training_data;
    delete prediction_data;

    // Create the retention time table:
    // RT : pointer to peptide sequence : rel. abundance
    for (UInt i = 0; i < peptidesVector.size(); i++)
    {
      PeptideSequences::const_iterator pit = sample_.getPeptideSequences().find( peptidesVector[i] );

      // check for rt outlier and rectify
      if (predicted_retention_times[i] < 0.0) predicted_retention_times[i] = 0.0;
      else if (predicted_retention_times[i] > 1.0) predicted_retention_times[i] = 1.0;

      rtTable.insert( make_pair( predicted_retention_times[i], pit ) );
    }

  #ifdef DEBUG_SIM
    cout << "---------------------------------------------------------" << endl;
    cout << "Content of retention time table:" << endl;
    for (RTTable::const_iterator rtTableRow = rtTable.begin();
          rtTableRow != rtTable.end();
          ++rtTableRow)
    {
      cout << rtTableRow->first << ", ";
      cout << ( (rtTableRow->second)->first ) << " , abundance: " << ( (rtTableRow->second)->second ) << endl;
    }
    cout << "---------------------------------------------------------" << endl;
  #endif
  }

  double LCMSSim::sampleModifications_(AASequence& aas, EmpiricalFormula& ef)
  {
    pair<ModTableIterator,ModTableIterator> range;

    for (Size i=0;i<aas.size();++i)
    {
      range = allowedMods_.equal_range( aas[i].getOneLetterCode() );

      for (ModTableIterator it = range.first; it != range.second; ++it)
      {
        double res = gsl_ran_flat(rand_gen_,0.0,1.0);

        if (res < it->second.abundance_)
        {
          // determine direction of shift (pos or negative)
          if (it->second.shift_)
            ef += it->second.formula_;
          else
            ef -= it->second.formula_;

          return it->second.abundance_;
        }

      }

    }

    return 0.0;
  }

  void LCMSSim::chooseElutionProfile_(ProductModel<2>& pm, const CoordinateType rt,const double scale)
  {
      Param p;
      double symmetry = gsl_ran_flat (rand_gen_, -100.0, 100.0);

      double width = gsl_ran_flat (rand_gen_, 5, 15);

      // Exponentially modified Gaussian
      const DoubleReal decay_stretch = 5.0;
      // it doesn't matter what bounding box I set here, it always cuts of the fronting !!! :-(
      p.setValue("bounding_box:min",rt - decay_stretch*(width+abs(symmetry)) );
      p.setValue("bounding_box:max",rt + decay_stretch*(width+abs(symmetry)) );
      p.setValue("interpolation_step", rt_sampling_ / 3.0);
      p.setValue("statistics:variance",1.0);
      p.setValue("statistics:mean",rt );
      p.setValue("emg:height",scale);
      p.setValue("emg:width",width);
      p.setValue("emg:symmetry",symmetry);
      p.setValue("emg:retention",rt);
      ElutionModel* elutionmodel = new ElutionModel();
      elutionmodel->setParameters(p);
      elutionmodel->setScalingFactor(scale);

      //----------------------------------------------------------------------

      // So far we had randomness in the parameters of the EMG.  Now let us add
      // some random ups and downs within the elution profile.

      // Hack away constness :-P   We know what we want.
      ElutionModel::ContainerType &data = const_cast<ElutionModel::ContainerType &>(elutionmodel->getInterpolation().getData());

      // Member distortion_ is the logarithm of the maximal factor which can be applied to
      // each point. But note that elution profile is smoothed again afterward.

      // first and last entry shall not be changed
      for ( UInt i = 1; i < data.size() - 1; ++i )
      {
        data[i] *= exp(gsl_ran_flat (rand_gen_, -distortion_, +distortion_));
      }
      // moving average filter (width 3), implemented very inefficiently (guess why!)
      if ( distortion_ != 0.0 ) // otherwise we want perfect EMG shape!
      {
        ElutionModel::ContainerType tmp;
        tmp.resize(data.size());
        for ( UInt i = 1; i < data.size() - 1; ++i )
        {
          tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
        }
        for ( UInt i = 1; i < data.size() - 1; ++i )
        {
          data[i] = tmp[i];
        }

        const int num_rounds = 10;
        for ( int rounds = 0; rounds < num_rounds; ++rounds)
        {
          data.swap(tmp);
          for ( UInt i = 1; i < data.size() - 1; ++i )
          {
            tmp[i] = ( data[i-1] + data[i] + data[i+1] ) / 3.0;
          }
          for ( UInt i = 1; i < data.size() - 1; ++i )
          {
            data[i] = tmp[i];
          }
        }

      }
      pm.setModel(0,elutionmodel);
  }

  void LCMSSim::addBaseline_()
  {
    double scale = param_.getValue("baseline_scaling");
    if (scale == 0.0) return;

    for (UInt i=0;i<exp_.size();++i)
    {
      for (UInt j=0;j<exp_[i].size();++j)
      {
        CoordinateType x = (exp_[i][j].getMZ() - minMapMZ_);
        if (x >= 1000.0) continue; // speed-up

        double b = gsl_ran_exponential_pdf(x,125.0);
        b *= scale;
        exp_[i][j].setIntensity( exp_[i][j].getIntensity() + b );
      }
    }
  }


  UInt LCMSSim::countBasicResidues_(const AASequence& seq) const
  {
    UInt count = 0;

    for (Size i = 0; i<seq.size(); ++i)
    {
      // check for basic residues
      // TODO: Currently we only consider "strong" basic residues. It might be a good idea to read this from a file.
      if (seq[i].getShortName() == "Arg" || seq[i].getShortName() == "Lys" || seq[i].getShortName() == "His")
      {
        ++count;
      }
    }

    return count;
  }


  void LCMSSim::samplePeptideModel_(const ProductModel<2>& pm,
                                                          const CoordinateType mz_start,  const CoordinateType mz_end,
                                                          CoordinateType rt_start,  CoordinateType rt_end )
  {
      // start and end points of the sampling are entirely arbitrary
      // and should be modified at some point

    // (cg) commented this out since it cuts off fronted elution profiles!!
    // (ost) Why should this happen ?
      if (rt_start <=0) rt_start = 0;

      LCMSmap::iterator exp_iter = exp_.RTBegin(rt_start);

      IntensityType intensity_sum = 0.0;
      vector< DPosition<2> > points;

      #ifdef DEBUG_SIM
      cout << "Sampling at: " << mz_start << " " << mz_end << " ";
      cout << rt_start << " " << rt_end << endl;
      #endif

      /// TODO: think of better error checking
      if(exp_iter == exp_.end() )
      {
        cout << "error ! " << endl; // ;-) should not happen
        return;
      }

      PointType point;

      Int start_scan = -5;
      Int end_scan  = -5;

      //UInt it = 0;
      //UInt pit = 0;

      for (CoordinateType rt = rt_start; rt < rt_end && exp_iter != exp_.end(); rt += rt_sampling_, ++exp_iter)
      {
        for (CoordinateType mz = mz_start; mz < mz_end; mz += msBinSize_)
        {
          //++it;

          point.setMZ(mz);
          point.setIntensity( pm.getIntensity( DPosition<2>( rt, mz) ) );

          if ( point.getIntensity() > 10.0)
          {
            if (start_scan == -5)
            {
              start_scan = exp_iter - exp_.begin();
              //cout << "start_scan: " << start_scan << endl;
            }

            if (! changed_scans_.at( exp_iter - exp_.begin() ) )
            {
              changed_scans_.at( exp_iter - exp_.begin() )  = true;
            }
            //++pit;

            // add m/z and itensity error (both Gaussian distributed)
            double it_err  = gsl_ran_gaussian(rand_gen_, (point.getIntensity() * intStdDevError_ ) ) + intMeanError_ ;

            // this is a quick fix only, should be improved to prevent simulation of negative intensities
            if (it_err < 0.0) it_err = fabs(it_err);

            point.setIntensity( point.getIntensity( ) + it_err );

            double mz_err = gsl_ran_gaussian(rand_gen_, mzStdDevError_) + mzMeanError_;
            point.setMZ( point.getMZ() + mz_err );

            intensity_sum += point.getIntensity();
            points.push_back( DPosition<2>( rt, mz) );		// store position
            exp_iter->push_back(point);
            //cout << "Sampling intensity: " << point.getIntensity() << endl;

            //update last scan
            end_scan = exp_iter - exp_.begin();
          }
        }
      }

      //cout << "End of sampling: " << it << " vs " << pit << endl;

      // do not set this here, because it might include low intensity points and is inconsistent with the convex hull
      //end_scan = exp_iter - exp_.begin();


      //cout << "end_scan: " << end_scan << endl;

      // This is a clear misuse of the Feature data structure
      // but at this point, we couldn't care less ;-)
      current_feature_.setQuality(0,start_scan);
      current_feature_.setQuality(1,end_scan);

      current_feature_.setIntensity(intensity_sum);
      // store convex hull
      current_feature_.getConvexHulls().clear();
      current_feature_.getConvexHulls().resize( current_feature_.getConvexHulls().size()+1);
      current_feature_.getConvexHulls()[ current_feature_.getConvexHulls().size()-1] = points;

      #ifdef DEBUG_SIM
      current_feature_.setModelDescription( ModelDescription<2>( &pm ) );
      #endif
  }

  void LCMSSim::addContaminants_()
  {
    double perc_contaminants = param_.getValue("perc_contaminants");
    double noise_it_mean       = param_.getValue("noise_int_mean");

    UInt num = (UInt) ceil(perc_contaminants * sample_.getPeptideSequences().size() / 100.0);
    if (num == 0) return;

    cout << "Adding contaminations." << endl;
    vector<EmpiricalFormula> vef;

    readFromContaminationFile_(vef);
		if (vef.size() == 0) return;

    ProductModel<2> pm;

    Param p1;
    Param p2;

    IsotopeModelGeneral* mz_model = 0;
    GaussModel* gauss2 = 0;

    mean_scaling_ /= ion_count_;
    mean_scaling_ *= 100.0;

    if (mean_scaling_ <= noise_it_mean)
    {
      mean_scaling_ = 2 * noise_it_mean;
    }

    cout << "Adding " << perc_contaminants << " % contaminants = "  << num << endl;
    cout << "Intensity scaling: " << mean_scaling_ << endl;

    for (UInt i=0;i<num;++i)
    {
      // Choose a contaminant at random
      UInt c_num = (UInt) gsl_ran_flat(rand_gen_, 0, vef.size());

      CoordinateType mass = vef[c_num].getMonoWeight();
      if (mass > maxMapMZ_ || mass < minMapMZ_) continue;

      // determine rt at random for non-peptdic contaminants (in a reasable range)
      CoordinateType rt = gsl_ran_flat(rand_gen_, 0, gradientTime_);
      rt = rt - 400.0 > 0 ? rt : (rt + abs(rt - 400.0) + 50.0);

      #ifdef DEBUG_SIM
      cout << "Inserting contaminant at mz " << mass << " and rt " << rt << endl;
      #endif

      p1.setValue("interpolation_step",0.0005);
      p1.setValue("statistics:mean",vef[c_num].getAverageWeight()); //  mass
      p1.setValue("isotope:stdev",peak_std_);
      p1.setValue("charge",1);

      mz_model = new IsotopeModelGeneral();
      mz_model->setParameters(p1);
      mz_model->setSamples(vef[c_num]);

      p2.setValue("bounding_box:min",rt-800);
      p2.setValue("bounding_box:max",rt+800);
      p2.setValue("statistics:variance",100.0);
      p2.setValue("statistics:mean",rt);
      gauss2 = new GaussModel();
      gauss2->setParameters(p2);
      gauss2->setSamples();

      pm.setModel(1,mz_model).setModel(0,gauss2);
      pm.setScale(mean_scaling_);

      current_feature_.setMZ(mass);
      current_feature_.setRT(rt);
      current_feature_.setCharge(1);

      samplePeptideModel_(pm, (mass-2.0),(mass+5.0), (rt-200),(rt+500) );

      if (current_feature_.getConvexHulls().begin()->getPoints().size() > 0)
      {
        current_feature_.setMetaValue("formula",vef[c_num].getString() );
        contaminants_.push_back(current_feature_);
      }
    }

  }

  void LCMSSim::removeDuplicatePoints_()
  {
#ifdef DEBUG_SIM
    cout << "Removing duplicates points in changed ....  " << flush;
#endif

    CoordinateType diff_mz = 0.0;
    PointType p;
    bool change = false;

    for(UInt i=0;i<exp_.size();++i)
    {
      // skip tiny or unchanged scans
      if (exp_[i].size() <= 50 || ! changed_scans_.at(i) ) continue;

      // copy Spectrum and clear Peaks
      LCMSmap::SpectrumType cont = exp_[i];
      cont.clear();

      exp_[i].sortByPosition();

      for (UInt j=0;j< (exp_[i].size() - 1);++j)
      {
          diff_mz = fabs(exp_[i][ (j+1) ].getMZ() - exp_[i][j].getMZ());

          if (diff_mz < msBinSize_)
          {
            change = true;
            // sum intensities
            CoordinateType it1 = exp_[i][ (j+1) ].getIntensity();
            CoordinateType it2 = exp_[i][ (j) ].getIntensity();
            CoordinateType it =  it1 + it2;
            p.setIntensity( it );

            // keep m/z of point with higher intensity
            CoordinateType mz1 = exp_[i][ (j+1) ].getMZ();
            CoordinateType mz2 = exp_[i][ (j) ].getMZ();
            CoordinateType mz = it1 > it2 ? mz1 : mz2;
            p.setMZ( mz );
            cont.push_back(p);

            ++j;
          }
          else
          {
            cont.push_back( exp_[i][j] );
          }
        }

      // reset flag
      changed_scans_.at(i) = false;

      // don't forget the last point
      if (!change) cont.push_back( exp_[i][ (exp_[i].size() - 1) ] );
      exp_[i] = cont;
    }

#ifdef DEBUG_SIM
    cout << "Done " << endl;
#endif
  }

  UInt LCMSSim::removeAllDuplicatePoints_()
  {
#ifdef DEBUG_SIM
    cout << "Removing all duplicate points....  " << flush;
#endif

    CoordinateType diff_mz = 0.0;
    PointType p;

    UInt count = 0;
    bool change = false;

    for(UInt i=0;i<exp_.size();++i)
    {
      exp_[i].sortByPosition();
      if (exp_[i].size() <= 2) continue;

      // copy Spectrum and remove Peaks ..
      LCMSmap::SpectrumType cont = exp_[i];
      cont.clear();

      for (UInt j=0;j< (exp_[i].size() - 1);++j)
      {
          diff_mz = fabs(exp_[i][ (j+1) ].getMZ() - exp_[i][j].getMZ());

          if (diff_mz < msBinSize_)
          {
            change = true;
            // sum intensities
            CoordinateType it1 = exp_[i][ (j+1) ].getIntensity();
            CoordinateType it2 = exp_[i][ (j) ].getIntensity();
            CoordinateType it =  it1 + it2;
            p.setIntensity( it );

            // keep m/z of point with higher intensity
            CoordinateType mz1 = exp_[i][ (j+1) ].getMZ();
            CoordinateType mz2 = exp_[i][ (j) ].getMZ();
            CoordinateType mz =  it1 > it2 ? mz1 : mz2;
            p.setMZ( mz );
            cont.push_back(p);

            ++j;
            ++count;
          }
          else
          {
            cont.push_back( exp_[i][j] );
          }
      }
      // don't forget the last one
      if (!change) cont.push_back( exp_[i][ (exp_[i].size() - 1) ] );
      exp_[i] = cont;
    }

#ifdef DEBUG_SIM
    cout << "Done " << endl;
    cout << "Count: " << count << endl;
#endif

    return count;
  }


  void LCMSSim::insertPeptideIon_(const EmpiricalFormula& ef, const CoordinateType rt, const ChargeType c, const double ab)
  {

    //TODO somewhere in here we could generate Tandem-MS scans:
    // what do we need for that?
    // - parent RT&mz
    // - AASequence (not given yet, needs to be added)
    // - ab(=Abundance)
    //
    // open Questions:	- where to store the resulting scan?! using exp_ is probably dangerous... introduce exp_tandem_? (unify them when storing mzData)
    //									- how is the resulting Tandem spectrum influenced by parent charge?
    //
    // TODO: define Interface
    // TODO: add parameters to this class to handle Tandem spectrum properties and labelling methods


    ProductModel<2> pm;
    Param p1;

    IntensityType scale = ab *1500; // was: 3000
    mean_scaling_ += scale;
    ++ion_count_;

    //TODO use H+ weight instead of 1*c ?
    CoordinateType mz = ( (ef.getMonoWeight() + c) / c) ;

    // TODO: this should not be necessary, check...
    allow_overlaps_ = (unsigned int) param_.getValue("allow_overlaps");
    if ( (allow_overlaps_ == 0) && checkForOverlaps_( DPosition<2>(mz,rt) ) )
    {
      return;
    }
    // don't show ions with m/z higher than the MS detection limit
    if (mz > maxMapMZ_ || mz < minMapMZ_) return;

#ifdef DEBUG_SIM
    cout << "Inserting ion with charge: " << c << " m/z " << mz << " rt " << rt << " rel. abundance " << ab;
    cout << " mass: " << ef.getMonoWeight() << endl;
#endif

    current_feature_.setMZ(mz);
    current_feature_.setRT(rt);
    current_feature_.setCharge(c);

    //TODO use H+ weight instead of 1*c ?
    p1.setValue("statistics:mean",((ef.getAverageWeight()+c)/c));
    p1.setValue("interpolation_step",0.001);
    p1.setValue("isotope:stdev",peak_std_);
    p1.setValue("charge",(int) c);

    IsotopeModelGeneral* isomodel = new IsotopeModelGeneral();
    isomodel->setSamples(ef);
    isomodel->setParameters(p1);

    chooseElutionProfile_(pm,rt,scale);
    pm.setModel(1,isomodel);
    pm.setScale(scale);

    // add peptide to global MS map
    // TODO: use flexible boundaries for rt/mz depending on abundance/charge?
    samplePeptideModel_(pm, (mz-2.5),(mz+5.0), (rt-160.0),(rt+280.0));
    if (current_feature_.getConvexHulls().begin()->getPoints().size() > 0)
    {
      features_.push_back(current_feature_);
    }

    // deletion of 1D models is done in destructor of productmodel, not necessary here
    //delete isomodel;

  }

  void LCMSSim::run()
  {
    cout << "Running simulation....  " << flush;

    // Predict retention times
    RTTable rtTable;
    predictRT_(rtTable);

    IonizationType MStype;

    UInt m = param_.getValue("ion_model");
    switch (m)
    {
      case 0: MStype = SIMPLE; break;
      case 1: MStype = ESI; break;
      default: MStype = ESI;
    }

    int num_scans = 1; // only one scan in case of 1D
    if(RTModelFile_ != "1D")
    {
      num_scans = (int) (gradientTime_ / rt_sampling_);
      num_scans += 400;	// add some noise scans
    }
    exp_.resize(num_scans);

    double scan_rt = rt_sampling_;
    for (UInt i=0; i<exp_.size(); ++i)
    {
      /// add some noise to retention times
      double n = gsl_ran_gaussian(rand_gen_, 0.05);
      exp_[i].setRT(scan_rt + n);
      scan_rt += rt_sampling_;
    }

    changed_scans_.resize(num_scans,false);

    CoordinateType rt  = 0.0;
    CoordinateType ab  = 0.0;

    // rt error
    CoordinateType rt_shift_mean  = param_.getValue("rt_shift_mean");
    CoordinateType rt_shift_stddev = param_.getValue("rt_shift_stddev");

    // Do we model modifications (yes/no) ?
    UInt ptms_on = (UInt)	param_.getValue("ptms_on");

    setLogType(CMD);
    startProgress(0, rtTable.size() , "simulation");
    UInt c = 0;

    // To convert AASequence to String
    stringstream str;
    // The peptde
    AASequence aas;
    // Relative ion abundances per charge state
    vector<double> charges;

    EmpiricalFormula ef;
    EmpiricalFormula ef_mod;

    for (RTTable::iterator cit = rtTable.begin();
          cit != rtTable.end();
          ++cit)
      {
        aas = ( (cit->second)->first );
        ab =  (cit->second)->second;								// peptide abundance
        rt = 400.0 + cit->first * gradientTime_; 		// rescale retention time

        // class AASequence does not have toString() method
        str << aas;
        // peptide sequence as string
        String aaseq_str = str.str();
        str.str("");

        setProgress(++c);

        if (c > 200 && (c % 20 == 0))
        {
          removeDuplicatePoints_();
        }

        switch (MStype)
        {
          case SIMPLE:
              // just one charge state (for debugging)
            charges.resize(2,1.0);
            break;

          case ESI:
            // we assume that the charge state for a peptide follows a Binomial distribution
            // TODO: think of better model
            UInt bas_c = countBasicResidues_(aas);
            charges.resize((bas_c+1),0);
            for (UInt i = 0; i<ab;++i)
            {
              cout << "ab : " << i << endl;
              cout << "max charge : " << bas_c << endl;
              //TODO: this is very dangerous! "i" has local and loop scope! (fixed)
              unsigned int bi = gsl_ran_binomial(rand_gen_,0.8,bas_c);
              ++charges[ bi ];
              cout << "Setting charge " << bi << " to " << charges[bi] << endl;
            }
            break;
        }

        // create LC-MS signal for each peptide ion
        // starting at c=1 (c=0 wouldnt show up in a spectrum)

        // rt error is independent of charge state
        // rt error is only modeled if we have a rt model
        CoordinateType rt_error = 0.0;
        if (RTModelFile_ != "none")
        {
           rt_error = gsl_ran_gaussian(rand_gen_, rt_shift_stddev) + rt_shift_mean;
        }
        for (UInt c=1;c<charges.size();++c)
        {
          if (charges[c] == 0) continue;

          // normalize relative ion abundances
          //charges[c] /= ab;
          ef = aas.getFormula();

          // retrieve abundance of modified peptide and create a second ion
          ef_mod = ef;
          double ab_mod = sampleModifications_(aas,ef_mod);

          if (ab_mod != 0.0 && ptms_on)
          {
#ifdef DEBUG_SIM
            cout << "Creating modified ion ! abundance : " << ab_mod;
            cout << "charge : " << c << endl;
#endif
            current_feature_.setMetaValue("sequence",aaseq_str);
            insertPeptideIon_(ef_mod,( rt+rt_error),c,(charges[c]*ab_mod));
            // set number of remaining "normal" instances
            charges[c] *= (1-ab_mod);
          }
#ifdef DEBUG_SIM
          cout << "Creating ion ! abundance : " << charges[c] << endl;
          cout << "charge : " << c << endl;
#endif
          current_feature_.setMetaValue("sequence",aaseq_str);
          insertPeptideIon_(ef,(rt+rt_error),c,charges[c]);
        }
        charges.clear();
      }

      endProgress();

      addContaminants_();
      addShotNoise_();

      UInt count = 0;
      do
      {
        count = removeAllDuplicatePoints_();
      } while (count != 0);

      addBaseline_();

      //cout << "Done." << endl;

  } // end of run()


  void LCMSSim::addShotNoise_()
  {
    // we model the amount of (background) noise as Poisson process
    // i.e. the number of noise data points per unit m/z interval follows a Poisson
    // distribution. Noise intensity is assumed to be Gaussian-distributed.

    double rate       = param_.getValue("noise_rate");
    double it_mean = param_.getValue("noise_int_mean");

    const UInt num_intervals = 100;
    CoordinateType interval_size = (maxMapMZ_ - minMapMZ_) / num_intervals;
    PointType point;

    cout << "Adding shot noise to spectra...." << endl;
    cout << "Interval size: "  << interval_size << " poisson rate: " << rate << endl;

    for (UInt i=0;i<exp_.size();++i)
    {

      for (UInt j=0;j<num_intervals;++j)
      {
        UInt counts = gsl_ran_poisson ( rand_gen_, rate);
        CoordinateType mz_lw = j * interval_size + minMapMZ_;
        CoordinateType mz_up = (j+1) * interval_size + minMapMZ_;

        for (UInt c=0; c<counts;++c)
        {
            CoordinateType mz  = gsl_ran_flat(rand_gen_, mz_lw, mz_up );
            CoordinateType it = gsl_ran_exponential(rand_gen_,it_mean);
            point.setIntensity(it);
            point.setMZ(mz);
            exp_[i].push_back(point);
        }

      }
    } // end of each scan

    exp_.updateRanges();

  }

  // Export spectrum as MzData file
  void LCMSSim::exportMzData(const String& filename)
  {
    MzDataFile().store(filename, exp_);
  }

  void LCMSSim::exportFeatureMap(const String& filename)
  {
    // write peptide lists

    Size i = filename.rfind(".");
    String xml_out = filename;
    String txt_out = filename;
    xml_out.replace(i,xml_out.size(),"_feature_list.featureXML");
    txt_out.replace(i,txt_out.size(),"_feature_list.txt");

    FeatureXMLFile().store(xml_out, features_);

    // text output
    ofstream out( txt_out.c_str() );
    // write header
    DateTime dt;
    String t; String d;
    dt.now();
    t = dt.getTime();
    d = dt.getDate();
    out << "# " << d << " " << t << endl;

    for (FeatureMap< >::const_iterator citer = features_.begin();
          citer != features_.end();
          ++citer)
    {
      out << citer->getPosition()[0] << " " << citer->getPosition()[1] << " " << citer->getIntensity();
      out << " " << citer->getCharge();
      out << " " << citer->getQuality(0);	// stores first scan
      out << " " << citer->getQuality(1);	// stores last scan
      out << " " << exp_[ (UInt) citer->getQuality(0) ].getRT();	// retrieve rt of first scan
      out << " " << exp_[ (UInt) citer->getQuality(1) ].getRT();	// retrieve rt of last scan
      out << " " << citer->getMetaValue("sequence");
      out << endl;
    }
    out.close();

    if (contaminants_.size() > 0)
    {
      // write metabolite lists
      xml_out = filename;
      txt_out = filename;
      xml_out.replace(i,xml_out.size(),"_contaminant_list.featureXML");
      txt_out.replace(i,txt_out.size(),"_contaminant_list.txt");

      FeatureXMLFile().store(xml_out, contaminants_);

      // text output
      out.open( txt_out.c_str() );
      // write header
      out << "# " << d << " " << t << endl;

      for (FeatureMap< >::const_iterator citer = contaminants_.begin();
            citer != contaminants_.end();
            ++citer)
      {
        out << citer->getPosition()[0] << " " << citer->getPosition()[1] << " " << citer->getIntensity();
        out << " " << citer->getCharge();
        out << " " << citer->getQuality(0);	// stores first scan
        out << " " << citer->getQuality(1);	// stores last scan
        out << " " << exp_[ (UInt) citer->getQuality(0) ].getRT();	// retrieve rt of first scan
        out << " " << exp_[ (UInt) citer->getQuality(1) ].getRT();	// retrieve rt of last scan
        out << " " << citer->getMetaValue("formula");
        out << endl;
      }
      out.close();
    }
  }

  void LCMSSim::updateMembers_()
  {
      gradientTime_ = param_.getValue("total_gradient_time");
      rt_sampling_ 	 = param_.getValue("rt_sampling");

      double tmp    =  param_.getValue("peak_fwhm");
      peak_std_     = (tmp / 2.355);			// Approximation for Gaussian-shaped signals
      msBinSize_      = param_.getValue("mz_sampling");

      maxMapMZ_    = param_.getValue("mz_upperlimit");
      minMapMZ_     = param_.getValue("mz_lowerlimit");

      mzMeanError_   = param_.getValue("mz_error_mean");
      mzStdDevError_ = param_.getValue("mz_error_stddev");

      intMeanError_   = param_.getValue("int_error_mean");
      intStdDevError_ = param_.getValue("int_error_stddev");

      String s             = param_.getValue("column_condition");

      allow_overlaps_ = (unsigned int) param_.getValue("allow_overlaps");

      if (s  == "poor")
      {
        distortion_ = 2.0;
        symmetry_down_ = -100;
        symmetry_up_ = +100;
      }
      else if (s == "medium")
      {
        distortion_ = 1.0;
        symmetry_down_ = -60;
        symmetry_up_ = +60;
      }
      else	// default is "good"
      {
        distortion_ = 0.0;
        symmetry_down_ = -15;
        symmetry_up_     = +15;
      }
      unsigned random_seed = (unsigned) param_.getValue("random_seed");

      #ifdef DEBUG_SIM
      std::cout << "random_seed: " << random_seed << std::endl;
      #endif

      if ( random_seed != 0 )
      {
        gsl_rng_set(rand_gen_,random_seed);
      }
  }

  bool LCMSSim::checkForOverlaps_(DPosition<2> pos)
  {
    // We do not really check for overlaps here, but only enforce
    // a certain minimum distance between features.......
    for (FeatureMap< >::const_iterator cit = features_.begin();
           cit != features_.end();
           ++cit)
    {
      double mz_dist = fabs(cit->getMZ() - pos[1]);
      double rt_dist = fabs(cit->getRT() - pos[0]);

      if (mz_dist < 6.0 && rt_dist < 350.0 )
        return true;

    }

    return false;
  }

  vector<String> LCMSSim::getValidColumnConditions()
  {
    vector<String> tmp;
    tmp.push_back("good");
    tmp.push_back("medium");
    tmp.push_back("poor");
    return tmp;
  }

  void LCMSSim::setRTModelFile(String RTModelFile)
  {
    RTModelFile_ = RTModelFile;
  }

  void LCMSSim::setSample(LCMSSample& s)
  {
    sample_ = s;
  }

} // namespace OpenMS
