// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Samuel Wein $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>



using namespace std;
using namespace OpenMS;


/** Type used to store mass and m/z values, currently double.
    Must be floating point. */
typedef double massT;

/** Type used to store intensity / abundance values, currently double.
    May be integer, but this will almost certainly require a lot of pain. */
typedef double intenT;


/** A histogram tracks the number of times a particular value occurs in a
    distribution.  Histograms are broken up into buckets; the bucket size
    is the width of the bucket.  Buckets are based at 0.*/
class Histogram : public RefBase {
  public:
    Histogram(double bucketSize) : _bucketSize(bucketSize), _inited(false), _numBuckets(0) {}
    Histogram(double bucketSize, double start, double end) : _bucketSize(bucketSize), _start(start),
        _inited(true), _numBuckets(0) { _buckets.resize(int(floor((end- start) / _bucketSize))); }
    // default copy constructor and assignment operator are okay.
    virtual ~Histogram() {}

    void add(double val);

    int getNumBuckets() const {
        return _numBuckets;
    }

    double getBucketStart(int i) const {
        return _start + (_bucketSize * double(i));
    }

    double getBucketEnd(int i) const {
        return _start + (_bucketSize * double(i+1));
    }

    int getBucketCount(int i) const {
        return _buckets[i];
    }

    ostream& write(ostream& out) const;
  private:
    double _bucketSize;
    double _start;
    bool _inited;
    vector<int> _buckets;
    countT _numBuckets;
};



struct HornBaseline {
    HornBaseline(MSSpectrum& spec, int numLeft, int numRight, int smoothOrder, int output) : _output(output) {
        savitskyGolayCoeffs(_savgol_coeffs, numLeft, numRight, 1, smoothOrder); //TODO IMPLEMENT IN OPENMS
        _min = spec.getMin();  //Minimum Mass
        massT max = spec.getMax(); // Maximum Mass
        massT baselineFrom = _min;
        int n = int(ceil((max - _min) / _baselinePeriod));
        for(int i = 0; i < n; ++i) {
            massT baselineTo = baselineFrom + _baselineRange;
            if (baselineTo > max) baselineTo = max;
            Baseln baseline = _detBaseline(spec.MZBegin(baselineFrom),spec.MZEnd(baselineTo)
            );
            _baseline.push_back(baseline);
//            if (_output) thrashOut << "Baseline " << baselineFrom <<" - " << baselineTo
//                << ": " << baseline.baselineInten << " ("
//                << (baseline.baselineInten / spec->getYMax() * 100.) << "%) " << baseline.noiseWidth
//                << endl;

            baselineFrom += _baselinePeriod;
        }
    }

    virtual ~HornBaseline() {}

    intenT getBaselineAt(MSSpectrum::Iterator& where) const {
        return _getBaseline((*where).first).baselineInten;
    }
    double getSignalToNoise(massT mass, intenT inten) const {
        Baseln baseline = _getBaseline(mass);
        return (inten - baseline.baselineInten) / baseline.noiseWidth;
    }
    virtual double getMaxSignalToNoiseIn(MSSpectrum& sub) const;
//    virtual ostream& write(ostream& thrashBaselineOut) const;
  private:

    struct Baseln {
        Baseln() : baselineInten(0), noiseWidth(0) {}
        Baseln(intenT _baselineInten, double _noiseWidth) : baselineInten(_baselineInten), noiseWidth(_noiseWidth) {}
        intenT baselineInten;
        double noiseWidth;
    };

    Baseln _getBaseline(massT where) const {

        int i = int(floor((where - _min) / _baselinePeriod));
        return _baseline[i];
    }


    Baseln _detBaseline(MSSpectrum::Iterator begin, MSSpectrum::Iterator end) { // DETERMINE BASELINE
        double maxint,minint;
        maxint=begin->getIntensity();
        minint=begin->getIntensity();
        for(MSSpectrum::Iterator i = begin; i != end; ++i) { //Iterate to find the max and min peaks
            if (i->getIntensity>maxint) maxint=i->getIntensity;
            if (i->getIntensity<minint) minint=i->getIntensity;
        }
        double bucketSize = ((maxint - minint) / _numBuckets);
        //cout << spec->getYMax() << " " << spec->getYMin() << " " << bucketSize << endl;
        if (bucketSize <= 0) return Baseln(0, 0);
        Histogram h(bucketSize);
        for(MSSpectrum::Iterator i = begin; i != end; ++i) { //CREATE histogram of each peak
            intenT inten = i->getIntensity();
            //if (inten < -10) throw InvalidArgumentException(String("Negative intensity of ")
            //    + inten + " at " + (*i).first + " (sample " + i.getIndex()
            //    + ") in " + spec->getName());
            h.add(inten);
        }

        int numBuckets = int(double(h.getNumBuckets()) * .25);
        int n = numBuckets + _savgol_coeffs.size();
        vector<double> d(n);
        for (int i = 0; i < numBuckets; ++i) {
            for (int j = 0; j <= i; ++j) {
                d[i] += h.getBucketCount(j);
            }
        }
        for (int i = numBuckets; i < n; ++i) d[i] = 0.; // zero pad.

        // find smoothed first derivative
        convolve(d, _savgol_coeffs, 1, d); //TODO implement

//        if (_output > 1) {
//            ofstream dout((String("out.baseline.") + spec->getXMin()).c_str());
//            dout << setprecision(8);
//            for(int i = 0; i < numBuckets; ++i) {
//                dout << (h.getBucketStart(i) + (h.getBucketEnd(i) - h.getBucketStart(i)) / 2)
//                    << " " << d[i] << endl;
//            }
//            dout.close();
//        }

        double maxden = 0.;
        double max = 0.;
        for(int i = 0; i < n; ++i) {
            if (d[i] > maxden) {
                maxden = d[i];
                max = h.getBucketStart(i) + (h.getBucketEnd(i) - h.getBucketStart(i)) / 2;
            }
        }
        double half = maxden / 2.;
        double starthalf = 0.;
        double endhalf = 0.;
        for (int i = 0; i < n; ++i) {
            if (d[i] > half) {
                starthalf = h.getBucketStart(i) + (h.getBucketEnd(i) - h.getBucketStart(i)) / 2;
                break;
            }
        }
        for (int i = n - 1; i >= 0; --i) {
            if (d[i] > half) {
                endhalf = h.getBucketStart(i) + (h.getBucketEnd(i) - h.getBucketStart(i)) / 2;
                break;
            }
        }
        return Baseln(max, endhalf - starthalf);
    }



    int _output;
    vector<double> _savgol_coeffs;

    massT _min;
    vector<Baseln> _baseline;

    static const int _numBuckets = 1000;
    static const massT _baselineRange;
    static const massT _baselinePeriod;

};




/**
    @page UTILS_THRASH THRASH

    @brief Deconvolute and deisotope an MzML file

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_THRASH.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_THRASH.html
*/
class TOPPTHRASH :
        public TOPPBase
{

public:
    TOPPTHRASH() :
        TOPPBase("THRASH", "Tool deconvolute and desisotope an MzML file.", false)
    {
    }

protected:
    void registerOptionsAndFlags_()
    {
        //registerStringOption_("in", "<sequence>","","Nucleic acid sequence",true,false);

        registerInputFile_("in_file","<infile>","","Input MzML file to THRASH\n",true,false);
        setValidFormats_("in_file", ListUtils::create<String>("MzML"));

        registerOutputFile_("out_file", "<outfile>", "", "The THRASHED spectra\n");
        setValidFormats_("out_file", ListUtils::create<String>("MzML"));
        registerDoubleList_("average_compositon","<avg_comp>",ListUtils::create<double>("9.75,12.25,3.75,7.0,1"),"The average composition of a monomer in order C,H,N,O,P");  // We default to RNA
        registerDoubleOption_("charge_unit","<charge_unit>",1.00727646,"The chargeunit, default is proton mass",false);
        registerIntOption_("mincharge","<mincharge>",1,"min charge to test",false);
        registerIntOption_("maxcharge","<maxcharge>",20,"max charge to test",false);
        registerDoubleOption_("signal_to_noise","<SNR>",4,"The minimum signal to noise ratio to identify a peak",false);
        registerDoubleOption_("RL_threshold","<rlThresh>",10,"The raw figure of merit threshold",false);
        //Resolution is gathered from MzML metadata
        registerFlag_("sub_mode","If present, subtract theoretical peak intensity as calculated from isotope distribution, otherwise subtract to baseline",true);
    }



    ExitCodes main_(int, const char**)
    {
        ProgressLogger progresslogger;
        progresslogger.setLogType(log_type_);

        String in_path(getStringOption_("in_file"));
        String out_path(getStringOption_("out_file"));
        int8_t maxCharge(getIntOption_("maxcharge"));
        int8_t minCharge(getIntOption_("mincharge"));
        DoubleList average_comp(getDoubleList_("average_composition"));
        double charge_unit(getDoubleOption_("charge_unit"));
        double snr(getDoubleOption_("signal_to_noise"));
        double rl_thresh(getDoubleOption_("RL_threshold"));
        bool sub_theor=getFlag_("sub_mode");



        // create MSExperiment
        MSExperiment<Peak1D> loaded_experiment;
        //for each spectra in the experiment

            //Find baseline
            //check that we are sorted, if not, sort
            //for each bucket
                //While we have a peak above threshold?
                    //extend window if we are near the end
                    //determine charge
                    //find possible isotope fits
                    //if we have a fit, better than FOM continue, otherwise iterate through all possible charge states
                    //If this isnt the same cluster we saw before add it as a cluster, otherwise finish with this bucket
                    //find the monoisotopic peak of the cluster and add it to our output spectrum

        //write output
        MzMLFile mtest;
        mtest.store(out_path, generated_exp);


        return EXECUTION_OK;
    }

    HornBaseline::HornBaseline(MSSpectrum &spec, int numLeft, int numRight, int smoothOrder, int output, ostream& thrashOut) : _output(output) { //What do we use spec for
        savitskyGolayCoeffs(_savgol_coeffs, numLeft, numRight, 1, smoothOrder); //TODO IMPLEMENT IN OPENMS
        _min = spec->getXMin();  //Minimum Mass
        massT max = spec->getXMax(); // Maximum Mass
        massT baselineFrom = _min;
        int n = int(ceil((max - _min) / _baselinePeriod));
        for(int i = 0; i < n; ++i) {
            massT baselineTo = baselineFrom + _baselineRange;
            if (baselineTo > max) baselineTo = max;
            Baseln baseline = _detBaseline(ref_cast<SpectrumConstRef>(
                spec->sub(baselineFrom, baselineTo)
            ));
            _baseline.push_back(baseline);
            if (_output) thrashOut << "Baseline " << baselineFrom <<" - " << baselineTo
                << ": " << baseline.baselineInten << " ("
                << (baseline.baselineInten / spec->getYMax() * 100.) << "%) " << baseline.noiseWidth
                << endl;

            baselineFrom += _baselinePeriod;
        }
    }



};

int main(int argc, const char** argv)
{
    TOPPTHRASH tool;
    return tool.main(argc, argv);
}
