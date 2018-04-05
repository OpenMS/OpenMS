// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_CONTAINER_H
#define OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_CONTAINER_H


#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <utility>
#include <functional>

#include <vector>
#include <set>
#include <map>


namespace OpenMS
{
  /**
        @ingroup Chemistry

        @brief Isotope distribution class

        Holds an isotope distribution with the weight value and according
        probability. Distribution can be add using the '+' or '+=' operators.

        The most important value which should be set is the max isotope value.
        This value can be set using the setMaxIsotope method. It is an upper
        bound for the number of isotopes which are calculated. E.g. if it is set
        to 3, only the first three isotopes, Monoisotopic mass, +1 and +2 are
        calculated.
        By default all possible isotopes are calculated, which leads to a large
        number of values, if the mass value is large!

        TODO indicate how this does not hold mass defect information, simply probabilities
        TODO indicate how this relates to other isotope models (Coarse ... )

    */
  class Element; 

  class OPENMS_DLLAPI IsotopeDistribution
  {
public:

    /// @name typedefs
    //@{
    /// container type, first holds the weight of the isotope, second the probability
    typedef Peak1D MassAbundance;
    typedef std::vector<MassAbundance> ContainerType;
    typedef ContainerType::iterator iterator;
    typedef ContainerType::iterator Iterator;
    typedef ContainerType::const_iterator const_iterator;
    typedef ContainerType::const_iterator ConstIterator;

    typedef ContainerType::reverse_iterator reverse_iterator;
    typedef ContainerType::reverse_iterator ReverseIterator;
    typedef ContainerType::const_reverse_iterator const_reverse_iterator;
    typedef ContainerType::const_reverse_iterator ConstReverseIterator;
    //@}


    enum Sorted { intensity, mass, undefined};


    /// @name Constructors and Destructors
    //@{
    /** Default constructor, note max_isotope must be set later
            @see setMaxIsotope(Size max_isotope)
    */
    IsotopeDistribution();

    /// Detailed constructor which sets the @p max_isotope
    explicit IsotopeDistribution(Size max_isotope);

    /// Copy constructor
    IsotopeDistribution(const IsotopeDistribution & isotope_distribution);

    /// Destructor
    virtual ~IsotopeDistribution();
    //@}

    /// @name Accessors
    //@{

    /// overwrites the container which holds the distribution using @p distribution
    void set(const ContainerType & distribution);

    /// returns the container which holds the distribution
    const ContainerType & getContainer() const;

    /// returns the maximal weight isotope which is stored in the distribution
    Peak1D::CoordinateType getMax() const;

    /// returns the minimal weight isotope which is stored in the distribution
    Peak1D::CoordinateType getMin() const;

    /// returns the size of the distribution which is the number of isotopes in the distribution
    Size size() const;

    /// clears the distribution and resets max isotope to 0
    void clear();

    /// remove intensities below the cutoff
    void trimIntensities(double cutoff);

    /// sort isotope distribution by intensity
    void sortByIntensity();

    /// sort isotope distribution by mass
    void sortByMass();

    /** @brief re-normalizes the sum of the probabilities of the isotopes to 1

            The re-normalisation is needed as in distributions with a lot of isotopes (and with high max isotope)
            the calculations tend to be inexact.
    */
    void renormalize();

    /** @brief Trims the right side of the isotope distribution to isotopes with a significant contribution.

            If the isotope distribution is calculated for large masses (and with high max isotope)
            it might happen that many entries contain only small numbers. This function can be
            used to remove these entries.

            Do consider normalising the distribution afterwards.
    */
    void trimRight(double cutoff);

    /** @brief Trims the left side of the isotope distribution to isotopes with a significant contribution.

            If the isotope distribution is calculated for large masses (and with high max isotope)
            it might happen that many entries contain only small numbers. This function can be
            used to remove these entries.

            Do consider normalising the distribution afterwards.
    */
    void trimLeft(double cutoff);


    
    bool isNormalized() const;


    double averageMass() const;

    bool isConvolutionUnit() const;
    //@}

    /// @name Operators
    //@{
    /// Assignment operator
    IsotopeDistribution & operator=(const IsotopeDistribution & isotope_distribution);

    /// equality operator, returns true if the @p isotope_distribution is identical to this, false else
    bool operator==(const IsotopeDistribution & isotope_distribution) const;

    /// inequality operator, returns true if the @p isotope_distribution differs from this, false else
    bool operator!=(const IsotopeDistribution & isotope_distribution) const;
    //@}

    /// @name Iterators
    //@{
    inline Iterator begin() { return distribution_.begin(); }

    inline Iterator end()   { return distribution_.end(); }

    inline ConstIterator begin() const { return distribution_.begin(); }

    inline ConstIterator end() const { return distribution_.end(); }

    inline ReverseIterator rbegin() { return distribution_.rbegin(); }

    inline ReverseIterator rend()   { return distribution_.rend(); }

    inline ConstReverseIterator rbegin() const { return distribution_.rbegin(); }

    inline ConstReverseIterator rend() const { return distribution_.rend(); }

    inline void insert(const Peak1D::CoordinateType& mass, const Peak1D::IntensityType& intensity)
    {
      distribution_.push_back(Peak1D(mass, intensity));
    }
    //@}

    /// @name Data Access Operators
    //@{
    /// operator which access a cell of the distribution and wraps it in SpectrumFragment struct
    Peak1D& operator[](const UInt64& index){ return distribution_[index];}

    //@}


protected:   

    /// sort wrapper of the distribution
    void sort_(std::function<bool(const MassAbundance& p1, const MassAbundance& p2)> sorter);
    /// takes a function as a parameter to transform the distribution
    void transform_(std::function<void(MassAbundance&)> lambda);

    /// stores the isotope distribution
    ContainerType distribution_;

    ///Holds the distribution sorted state
    Sorted sort_type;
  };


} // namespace OpenMS

#endif // OPENMS_CHEMISTRY_ISOTOPEDISTRIBUTION_CONTAINER_H
