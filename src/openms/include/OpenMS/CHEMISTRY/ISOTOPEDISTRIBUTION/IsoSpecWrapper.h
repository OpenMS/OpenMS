// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/KERNEL/Peak1D.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>


// Override IsoSpec's use of mmap whenever it is available
#define ISOSPEC_GOT_SYSTEM_MMAN false
#define ISOSPEC_GOT_MMAN false

#include <OpenMS/../../thirdparty/IsoSpec/IsoSpec/isoSpec++.h>


namespace OpenMS
{

  /**
    * @brief Interface to the IsoSpec algorithm.
    * 
    * Provides an interface to the IsoSpec algorithm.
    *
    * @code
    * Łącki MK, Startek M, Valkenborg D, Gambin A.
    * IsoSpec: Hyperfast Fine Structure Calculator.
    * Anal Chem. 2017 Mar 21;89(6):3272-3277. doi: 10.1021/acs.analchem.6b01459.
    * @endcode
    *
    **/
  class OPENMS_DLLAPI IsoSpecWrapper
  {
  /**
   * @brief Interface for the IsoSpec algorithm - a generator of infinitely-resolved theoretical spectra.
   */

public:

    /**
      * @brief Run the algorithm
      *
      * This method will run the algorithm with parameters as set up by the constructor. It will return an
      * IsotopeDistribution containing the observed configurations. The configurations are explicitly stored
      * in memory, which may become a problem when considering some especially large distributions. If this,
      * or (a rather small) performance overhead is a concern, then the generator methods (below) should be
      * used instead.
      *
      * This method is provided for convience. As calling that method invalidates the object (the method should
      * not be called again, nor anythong other than destroying the object should be done with it), the most common
      * usage pattern of IsoSpecWrapper classes with the run method is:
      *
      * IsotopeDistribution dist = IsoSpecSubclass(...).run();
      * do something with dist;
      *
      * @note Calling this method invalidates the object! In future versions this limitation will be removed.
      *
      * @note This method should not be mixed with the generator methods - a given object should only ever have
      * its run method used, or only its generator methods used.
      *
      **/
    virtual IsotopeDistribution run() = 0;

    /**
     * @brief Move the generator to a next isotopologue
     *
     * Advance the internal generator to the next isotopologue. The value returned determines whether the
     * generator has been exhausted. It is invalid to call any other generator methods before the first call
     * to nextConf(), as well as after it returns false.
     *
     * A common, correct usage pattern would be:
     *
     * IsotopeDistributionSubclass isotopeDist(...);
     *
     * while(isotopeDist.nextConf())
     * {
     *     Peak1D conf = isotopeDist.getConf(); // and/or GetMass, GetProb, etc.
     *     do some computations on that conf;
     * }
     *
     * @returns A boolean value stating whether the generator has been exhausted.
     */
    virtual bool nextConf() = 0;

    /**
     * @brief Obtain the current isotopologue
     * @return The current isotopologue as a Peak1D
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual Peak1D getConf() = 0;

    /**
     * @brief Obtain the mass of the current isotopologue
     * @return The mass of the current isotopologue
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getMass() = 0;

    /**
     * @brief Obtain the intensity (probability, relative peak height) of the current configuration
     * @return The intensity (probability) of the current isotopologue
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getIntensity() = 0;

    /**
     * @brief Obtain the natural logarithm of the intensity (probability, relative peak height) of the current configuration
     * @return The natural logarithm of intensity (probability) of the current isotopologue
     * @note It is invalid (undefined results) to call this method before the first call to nextConf(), or after it returns false
     */
    virtual double getLogIntensity() = 0;

    /**
     * @brief Destructor
     */
    virtual inline ~IsoSpecWrapper() {};
  };

  class OPENMS_DLLAPI IsoSpecThresholdWrapper : public IsoSpecWrapper
  {
    /**
     * @brief Provides a threshold-based generator of isotopologues: generates all isotopologues
     *        more probable than a given probability threshold.
     *
     * This is the simplest generator - most users will however want to use IsoSpecTotalProbWrapper.
     * The reason for it is that when thresholding by peak intensity one has no idea how far the obtained
     * spectrum is from a real spectrum. For example, consider human insulin: if the threshold is set at
     * 0.1 of the most intense peak, then the total amount covered will be 99.7% for water, 82% for substance P,
     * only 60% for human insulin, and 23% for titin.
     * For a threshold of 0.01, the numbers will be: still 99.7% for water, 96% for substance P, 88% for human insulin
     * and 72% for titin (it also took 5 minutes on an average notebook computer to process the 17 billion configurations
     * involved).
     *
     * As you can see the threshold does not have a straightforward correlation to the accuracy of the final spectrum
     * obtained - and accuracy of final spectrum is often what the user is interested in. The IsoSpecProbabilityWrapper
     * provides a way to directly parametrise based on the desired accuracy of the final spectrum - and should be used
     * instead in most cases. The tradeoff is that it's (slightly) slower than Threshold algorithm. This speed gap will
     * be dramatically improved with IsoSpec 2.0.
     *
     * @note The isotopologues are NOT guaranteed to be generated in any particular order.
     */

public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param threshold Intensity threshold: will only compute peaks above this threshold
      * @param absolute Whether the threshold is absolute or relative (relative to the most intense peak)
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
  IsoSpecThresholdWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double threshold,
             bool absolute);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
  IsoSpecThresholdWrapper(const EmpiricalFormula& formula, double threshold, bool absolute);

  virtual IsotopeDistribution run() override final;

  virtual inline bool nextConf() override final { return ITG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(ITG.mass(), ITG.prob()); };
  virtual inline double getMass() override final { return ITG.mass(); };
  virtual inline double getIntensity() override final { return ITG.prob(); };
  virtual inline double getLogIntensity() override final { return ITG.lprob(); };


protected:
  IsoSpec::IsoThresholdGenerator ITG;
  };




  class OPENMS_DLLAPI IsoSpecTotalProbWrapper : public IsoSpecWrapper
  {
  /**
   * @brief Generate an optimal p-set of configurations for a given p (that is, a spectrum which smallest
   *        amount of isotopologues that has p accuracy)
   *
   * An optimal p-set of isotopologues is the smallest set of isotopologues that, taken together, cover at
   * least p of the probability space (that is, their probabilities sum up to at least p). This means that
   * the computed spectrum is accurate to at least degree p, and that the L1 distance between the computed
   * spectrum and the true spectrum is less than 1-p.
   *
   * This is the method most users will want: the p parameter directly controls the accuracy of results.
   *
   * @note The configurations are not guaranteed to be returned in any particular order.
   */
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      * @param total_prob Total coverage of probability space desired (e.g. 0.99)
      *
      * @note This constructor is only useful if you need to define non-standard abundances
      *       of isotopes, for other uses the one accepting EmpiricalFormula is easier to use.
      *
      **/
  IsoSpecTotalProbWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities,
             double total_prob);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
  IsoSpecTotalProbWrapper(const EmpiricalFormula& formula, double total_prob);

  virtual IsotopeDistribution run() override final;

  virtual inline bool nextConf() override final { return ILG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(ILG.mass(), ILG.prob()); };
  virtual inline double getMass() override final { return ILG.mass(); };
  virtual inline double getIntensity() override final { return ILG.prob(); };
  virtual inline double getLogIntensity() override final { return ILG.lprob(); };

protected:
  IsoSpec::IsoLayeredGenerator ILG;
  };



  class OPENMS_DLLAPI IsoSpecOrderedGeneratorWrapper : public IsoSpecWrapper
  {
  /**
   * @brief Generate the stream of configurations, ordered from most likely to least likely.
   *
   * This generator walks through the entire set of isotopologues of a given molecule, allowing
   * the user to terminate the search on the fly. The returned isotopologues are guaranteed to
   * be generated in order of descending probability (unlike the previous two generators which
   * make no such guarantees).
   *
   * This causes the algorithm to run in O(n*log(n)) and means that is it much slower than the
   * previous two.
   *
   * You should only use this generator if you don't know up-front when to stop the walk through
   * the configuration space, and need to make the decision on the fly. If you know the threshold
   * or the total probability needed, and only need the configurations sorted, it will be much
   * faster to generate them using one of the previous algorithms and sort them afterwards.
   */
public:
    /**
      * @brief Constructor
      *
      * @param isotopeNumbers A vector of how many isotopes each element has, e.g. [2, 2, 3])
      * @param atomCounts How many atoms of each we have [e.g. 12, 6, 6 for Glucose]
      * @param isotopeMasses Array with the individual elements isotopic masses
      * @param isotopeProbabilities Array with the individual elements isotopic probabilities
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const std::vector<int>& isotopeNumbers,
             const std::vector<int>& atomCounts,
             const std::vector<std::vector<double> >& isotopeMasses,
             const std::vector<std::vector<double> >& isotopeProbabilities);

    /**
      * @brief Setup the algorithm to run on an EmpiricalFormula
      *
      **/
  IsoSpecOrderedGeneratorWrapper(const EmpiricalFormula& formula);

  virtual IsotopeDistribution run() override final
  { throw std::logic_error("There is no stop condition in OrderedGenerator - therefore it only makes sense to use it as a generator"); };

  virtual inline bool nextConf() override final { return IOG.advanceToNextConfiguration(); };
  virtual inline Peak1D getConf() override final { return Peak1D(IOG.mass(), IOG.prob()); };
  virtual inline double getMass() override final { return IOG.mass(); };
  virtual inline double getIntensity() override final { return IOG.prob(); };
  virtual inline double getLogIntensity() override final { return IOG.lprob(); };

protected:
  IsoSpec::IsoOrderedGenerator IOG;
  };

}

