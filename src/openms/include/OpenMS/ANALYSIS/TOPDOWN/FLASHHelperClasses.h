// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#ifdef _MSC_VER
using uint = unsigned int; // POSIX uint not available on MSVC; was provided transitively via Qt
#endif
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <boost/dynamic_bitset.hpp>
#include <fstream>
#include <limits>
#include <unordered_set>

namespace OpenMS
{
  class PeakGroup;
  class ResidueModification;
  class SpectralDeconvolution;

  /**
   * @brief Wrapper struct for all the structs needed by the FLASHDeconv
   * The following structures/classes are defined: PrecalculatedAveragine, MassFeature, IsobaricQuantities, LogMzPeak,
   * Tag, and DAG for tagging.
   * i) PrecalculatedAveragine - to match observed isotopic envelope against theoretical one, theoretical envelope from
   * averagine model should be quickly calculated. To do so, precalculate averagines for different masses at the beginning of FLASHDeconv runs
   * ii) MassFeature - the feature of the deconvolved masses (instead of m/zs for MS features).
   * iii) IsobaricQuantities - the quantities from isobaric quantification (implemented within FLASHDeconv)
   * iv) LogMzPeak - Log transformed peak from original peak. Contains information such as charge, isotope index, and
   * uncharged mass.
   * v) Tag - the sequence tag generated on deconvolved spectra in FLASHTaggerAlgorithm
   * vi) DAG - the graph for sequence tagging in FLASHTaggerAlgorithm and for proteoform characterization in FLASHExtenderAlgorithm
   * @see SpectralDeconvolution
   */

  struct OPENMS_DLLAPI FLASHHelperClasses
  {
    /// @brief Averagine patterns pre-calculated for speed up. Other variables are also calculated for fast cosine calculation
    class OPENMS_DLLAPI PrecalculatedAveragine
    {
    private:
      /// isotope distributions for different (binned) masses
      std::vector<IsotopeDistribution> isotopes_;
      /// L2 norms_ for masses - for quick isotope cosine calculation
      std::vector<double> norms_;
      /// mass differences between average mass and monoisotopic mass
      std::vector<double> average_mono_mass_difference_;
      /// mass differences between most abundant mass and monoisotopic mass
      std::vector<double> abundant_mono_mass_difference_;

      std::vector<double> snr_mul_factor_;
      /// Isotope start indices: isotopes of the indices less than them have very low intensities
      std::vector<int> left_count_from_apex_;
      /// Isotope end indices: isotopes of the indices larger than them have very low intensities
      std::vector<int> right_count_from_apex_;
      /// most abundant isotope index
      std::vector<Size> apex_index_;

      /// max isotope index
      Size max_isotope_index_;
      /// mass interval for calculation
      double mass_interval_;
      /// min mass for calculation
      double min_mass_;
      /// calculate the mass bin index from mass
      Size massToIndex_(double mass) const;

      friend class OpenMS::PeakGroup;
      friend class OpenMS::SpectralDeconvolution;

      /// borrow the cached distribution for read-only scoring (same bin as get(); valid while this model lives)
      const IsotopeDistribution& getReference_(double mass) const
      {
        return isotopes_[massToIndex_(mass)];
      }

    public:
      /// default constructor
      PrecalculatedAveragine() = default;

      /**
       @brief constructor with parameters such as mass ranges and bin size.
       @param min_mass the averagine distributions will be calculated from this min_mass
       @param max_mass to the max_mass
       @param delta with the bin size delta
       @param generator this generates (calculates) the distributions
       @param use_RNA_averagine if set, nucleotide-based isotope patters are calculated
       @param decoy_iso_distance if set to a positive value, nonsensical isotope patterns are generated - the distance between isotope = decoy_iso_distance * normal distance.
    */
      PrecalculatedAveragine(double min_mass, double max_mass, double delta, CoarseIsotopePatternGenerator& generator, bool use_RNA_averagine, double decoy_iso_distance = -1);

      /// copy constructor
      PrecalculatedAveragine(const PrecalculatedAveragine&) = default;

      /// move constructor
      PrecalculatedAveragine(PrecalculatedAveragine&& other) noexcept = default;

      /// copy assignment operator
      PrecalculatedAveragine& operator=(const PrecalculatedAveragine& pc) = default;

      /// move assignment operator
      PrecalculatedAveragine& operator=(PrecalculatedAveragine&& pc) noexcept = default;

      /// destructor
      ~PrecalculatedAveragine() = default;


      /// get distribution for input mass. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      IsotopeDistribution get(double mass) const;

      /// get max isotope index
      size_t getMaxIsotopeIndex() const;

      /// set max isotope index
      void setMaxIsotopeIndex(int index);

      /// get isotope distance (from apex to the left direction) to consider. This is specified in the constructor. index. If input mass exceeds the
      /// maximum mass (specified in constructor), output for the maximum mass
      Size getLeftCountFromApex(double mass) const;

      /// get isotope distance (from apex to the right direction) to consider. This is specified in the constructor. index. If input mass exceeds the
      /// maximum mass (specified in constructor), output for the maximum mass
      Size getRightCountFromApex(double mass) const;

      /// get index of most abundant isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getApexIndex(double mass) const;

      /// get index of last isotope. If input mass exceeds the maximum mass (specified in constructor), output for the maximum mass
      Size getLastIndex(double mass) const;

      /// get mass difference between avg and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for the maximum
      /// mass
      double getAverageMassDelta(double mass) const;

      /// get mass difference between most abundant mass and mono masses. If input mass exceeds the maximum mass (specified in constructor), output for
      /// the maximum mass
      double getMostAbundantMassDelta(double mass) const;

      /// get the SNR multiplication factor - used for quick SNR calculation
      double getSNRMultiplicationFactor(double mass) const;
    };

    /// Mass feature (Deconvolved masses in spectra are traced by Mass tracing to generate mass features - like LC-MS features).
    struct OPENMS_DLLAPI MassFeature
    {
    public:
      /// feature index;
      uint index;
      /// the trace calculated from the masses
      MassTrace mt;
      /// per charge and isotope intensities
      std::vector<float> per_charge_intensity;
      std::vector<float> per_isotope_intensity;
      /// isotope offset between deconvolved masses and mass feature
      int iso_offset;
      /// representative scan number, min and max scan numbers, and representative charge
      int scan_number, min_scan_number, max_scan_number, rep_charge;
      /// average mass
      double avg_mass;
      /// minimum maximum charges, and number of charges
      int min_charge, max_charge, charge_count;
      double isotope_score, qscore;
      double rep_mz;
      bool is_decoy;
      uint ms_level;
      /// features are compared
      bool operator<(const MassFeature& a) const
      {
        return avg_mass < a.avg_mass;
      }
      bool operator>(const MassFeature& a) const
      {
        return avg_mass > a.avg_mass;
      }
      bool operator==(const MassFeature& other) const
      {
        return avg_mass == other.avg_mass;
      }
    };

    /// Isobaric quantities.
    struct OPENMS_DLLAPI IsobaricQuantities
    {
    public:
      int scan;
      double rt;
      double precursor_mz;
      double precursor_mass;
      std::vector<double> quantities;
      std::vector<double> merged_quantities;
      bool empty() const;
    };

    /// log transformed peak. After deconvolution, all necessary information from deconvolution such as charge and isotope index is stored.
    class OPENMS_DLLAPI LogMzPeak
    {
    public:
      /// original peak mz
      double mz = 0;
      /// original peak intensity
      float intensity = 0;
      /// log transformed mz
      double logMz = -1000;
      /// determined mass after deconvolution. NOT monoisotopic but only decharged
      double mass = .0;
      /// absolute charge (in case negative, is_positive is set to false)
      int abs_charge = 0;
      /// is positive mode
      bool is_positive = true;
      /// isotope index
      int isotopeIndex = -1;

      /// default constructor
      LogMzPeak() = default;

      /**
        @brief constructor from Peak1D.
        @param peak the original spectral peak
        @param positive determines the charge carrier mass. Can be obtained by getChargeMass(true) for positive mode (Constants::PROTON_MASS_U) and
        getChargeMass(false) for negative mode
        (-Constants::PROTON_MASS_U)
      */
      explicit LogMzPeak(const Peak1D& peak, bool positive);

      /// constructor from Peak1D with its log m/z already computed (getLogMz(peak.getMZ(), positive))
      LogMzPeak(const Peak1D& peak, bool positive, double log_mz) :
          mz(peak.getMZ()), intensity(peak.getIntensity()), logMz(log_mz), abs_charge(0), is_positive(positive), isotopeIndex(0)
      {
      }

      /// copy constructor
      LogMzPeak(const LogMzPeak&) = default;

      /// destructor
      ~LogMzPeak() = default;

      /// get uncharged mass of this peak. It is NOT a monoisotopic mass of a PeakGroup, rather the mass of each LogMzPeak. Returns 0 if no
      /// charge is set
      double getUnchargedMass() const;

      /// log mz values are compared
      bool operator<(const LogMzPeak& a) const;
      bool operator>(const LogMzPeak& a) const;
      bool operator==(const LogMzPeak& other) const;
    };

    class OPENMS_DLLAPI MassTag
    {
    public:
      /// constructor
      explicit MassTag(const std::string& full_seq, const std::vector<int>& positions, int scan);

      /// copy constructor
      MassTag(const MassTag&) = default;

      /// destructor
      ~MassTag() = default;

      bool operator<(const MassTag& a) const;
      bool operator>(const MassTag& a) const;
      bool operator==(const MassTag& other) const;

      [[nodiscard]] const std::string& getSequence() const;
      [[nodiscard]] double getNtermMass() const;
      [[nodiscard]] double getCtermMass() const;
      Size getLength() const;
      int getScore() const;
      int getScore(int pos) const;
      int getScan() const;
      int getIndex() const
      {
        return index_;
      }

      void setIndex(int i)
      {
        index_ = i;
      }

      float getRetentionTime() const
      {
        return rt_;
      }

      void setRetentionTime(float i)
      {
        rt_ = i;
      }

      std::string toString() const;
      const std::vector<double>& getMzs() const;

    private:
      double n_mass_ = -1, c_mass_ = -1;
      std::string seq_;
      float rt_ = -1;
      int scan_;
      std::vector<int> positions_;
      std::vector<double> masses_;
      int index_;
      Size length_;
    };


    /// Sequence tag used for FLASHTaggerAlgorithm. No mass gap is allowed in the seq.
    /// The mass gap containing tag should be enumerated into multiple Tag instances from outside.
    class OPENMS_DLLAPI Tag
    {
    public:
      /// constructor
      explicit Tag(const std::string& seq, double n_mass, double c_mass, const std::vector<double>& mzs, const std::vector<int>& scores, int scan);

      /// constructor of a gapped tag: @p mzs are the masses at the gap boundaries only, so the residue length is given explicitly
      explicit Tag(std::string seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores, int scan, Size length);

      /// copy constructor
      Tag(const Tag&) = default;

      /// move constructor (the declared copy constructor would otherwise make every move a copy)
      Tag(Tag&&) noexcept = default;

      Tag& operator=(const Tag&) = default;
      Tag& operator=(Tag&&) noexcept = default;

      /// destructor
      ~Tag() = default;

      bool operator<(const Tag& a) const;
      bool operator>(const Tag& a) const;
      bool operator==(const Tag& other) const;

      [[nodiscard]] const std::string& getSequence() const;
      [[nodiscard]] const std::string& getUppercaseSequence() const;
      [[nodiscard]] double getNtermMass() const;
      [[nodiscard]] double getCtermMass() const;
      /// number of residues the tag spans (for a gapped tag more than the number of mass gaps)
      Size getLength() const;
      /// true if consecutive masses of the tag may be more than one residue apart
      bool isGapped() const
      {
        return length_ + 1 != mzs_.size();
      }
      int getScore() const;
      int getScore(int pos) const;
      int getScan() const;
      int getIndex() const
      {
        return index_;
      }

      void setIndex(int i)
      {
        index_ = i;
      }

      float getRetentionTime() const
      {
        return rt_;
      }

      void setRetentionTime(float i)
      {
        rt_ = i;
      }

      std::string toString() const;
      const std::vector<double>& getMzs() const;

    private:
      std::string seq_, upper_seq_;
      double n_mass_ = -1, c_mass_ = -1;
      float rt_ = -1;
      std::vector<double> mzs_;
      std::vector<int> scores_;
      int scan_;
      int index_;
      Size length_;
    };

    /// A generic Directed Acyclic Graph structure for FLASHTaggerAlgorithm and FLASHExtensionAlgorithm.
    /// It is a 1-D graph so multiple dimension graph should be changed into 1-D to use this structure.
    /**
       @brief A set of vertex indices for graphs whose vertex space is huge but sparsely visited. A paged bitset: pages of 2^15 vertices
       (4 KB) are created when first touched, so memory follows the touched part of the space, and a lookup is a page lookup (the last
       page is cached; vertex indices of one search are clustered) plus a bit test. Pages are kept in blocks of 16 that never move (one
       growing array copied all pages at every reallocation); clear() keeps the blocks for the next use, and a page is zeroed when it is
       handed out again (a new block is zero already).
    */
    class OPENMS_DLLAPI VisitedVertexSet
    {
    public:
      /// true if vertex was not in the set
      bool insert(Size vertex)
      {
        uint64_t* page = page_(vertex >> page_bits_, true);
        const Size bit = vertex & (page_size_ - 1);
        uint64_t& word = page[bit >> 6];
        const uint64_t mask = uint64_t(1) << (bit & 63);
        if (word & mask) return false;
        word |= mask;
        ++count_;
        return true;
      }

      bool contains(Size vertex) const
      {
        const uint64_t* page = const_cast<VisitedVertexSet*>(this)->page_(vertex >> page_bits_, false);
        if (page == nullptr) return false;
        const Size bit = vertex & (page_size_ - 1);
        return (page[bit >> 6] >> (bit & 63)) & 1;
      }

      void clear()
      {
        std::fill(table_keys_.begin(), table_keys_.end(), empty_);
        page_count_ = 0;
        count_ = 0;
        last_page_id_ = empty_;
      }

      Size size() const
      {
        return count_;
      }

    private:
      static constexpr Size page_bits_ = 15;
      static constexpr Size page_size_ = Size(1) << page_bits_;
      static constexpr Size page_words_ = page_size_ / 64;
      static constexpr Size empty_ = std::numeric_limits<Size>::max();
      static constexpr Size block_bits_ = 4;
      static constexpr Size block_pages_ = Size(1) << block_bits_;

      std::vector<std::vector<uint64_t>> blocks_; // pages, block_pages_ per block, page_words_ words each
      Size used_pages_ = 0;                        // pages handed out since the blocks were allocated (those below are zeroed on reuse)
      std::vector<Size> table_keys_;    // page id per slot (open addressing)
      std::vector<Size> table_values_;  // page number per slot
      Size table_mask_ = 0;
      Size page_count_ = 0;
      Size count_ = 0;
      Size last_page_id_ = empty_;
      Size last_page_number_ = 0;

      static Size slot_(Size id, Size mask)
      {
        return (Size)((id * 0x9E3779B97F4A7C15ULL) >> 24) & mask;
      }

      uint64_t* pageData_(Size number)
      {
        return blocks_[number >> block_bits_].data() + (number & (block_pages_ - 1)) * page_words_;
      }

      uint64_t* page_(Size id, bool create)
      {
        if (id == last_page_id_) return pageData_(last_page_number_);
        if (! table_keys_.empty())
        {
          for (Size i = slot_(id, table_mask_);; i = (i + 1) & table_mask_)
          {
            if (table_keys_[i] == id)
            {
              last_page_id_ = id;
              last_page_number_ = table_values_[i];
              return pageData_(last_page_number_);
            }
            if (table_keys_[i] == empty_) break;
          }
        }
        if (! create) return nullptr;
        if ((page_count_ + 1) * 2 > table_keys_.size()) grow_table_();
        Size i = slot_(id, table_mask_);
        while (table_keys_[i] != empty_) i = (i + 1) & table_mask_;
        table_keys_[i] = id;
        table_values_[i] = page_count_;
        if ((page_count_ >> block_bits_) == blocks_.size()) blocks_.emplace_back(block_pages_ * page_words_, uint64_t(0));
        uint64_t* page = pageData_(page_count_);
        if (page_count_ < used_pages_) std::fill(page, page + page_words_, uint64_t(0));
        else used_pages_ = page_count_ + 1;
        last_page_id_ = id;
        last_page_number_ = page_count_++;
        return page;
      }

      void grow_table_()
      {
        std::vector<Size> old_keys, old_values;
        old_keys.swap(table_keys_);
        old_values.swap(table_values_);
        const Size n = old_keys.empty() ? 64 : old_keys.size() * 2;
        table_keys_.assign(n, empty_);
        table_values_.assign(n, 0);
        table_mask_ = n - 1;
        for (Size k = 0; k < old_keys.size(); k++)
        {
          if (old_keys[k] == empty_) continue;
          Size i = slot_(old_keys[k], table_mask_);
          while (table_keys_[i] != empty_) i = (i + 1) & table_mask_;
          table_keys_[i] = old_keys[k];
          table_values_[i] = old_values[k];
        }
      }
    };

    class OPENMS_DLLAPI DAG
    {
    public:
      DAG()
      {
      }
      explicit DAG(Size vertice_count): vertex_count_(vertice_count)
      {
      }

      Size size() const
      {
        return vertex_count_;
      }

      /// same as the bitset version, for callers whose vertex space is huge but sparsely visited (FLASHExtender)
      bool addEdge(Size vertex1, Size vertex2, VisitedVertexSet& visited)
      {
        if (vertex1 >= vertex_count_ || vertex2 >= vertex_count_) return false;
        if (! visited.contains(vertex2)) return false;
        visited.insert(vertex1);
        edges_.emplace_back(vertex1, vertex2);
        return true;
      }

      /// addEdge for a @p vertex2 the caller knows to be in @p visited: the same without that lookup (FLASHExtender adds the edges of a
      /// search from its start vertex, which it tests once)
      bool addEdgeFromVisited(Size vertex1, Size vertex2, VisitedVertexSet& visited)
      {
        if (vertex1 >= vertex_count_ || vertex2 >= vertex_count_) return false;
        visited.insert(vertex1);
        edges_.emplace_back(vertex1, vertex2);
        return true;
      }

      bool addEdge(Size vertex1, Size vertex2, boost::dynamic_bitset<>& visited)
      {
        if (vertex1 >= vertex_count_ || vertex2 >= vertex_count_) return false;
        if (! visited[vertex2]) return false;
        visited[vertex1] = true;
        edges_.emplace_back(vertex1, vertex2); //
        return true;
      }

      /**
       * Simple path finding algorithm on DAG
       * @param source source
       * @param sink sink
       * @param all_paths the resulting paths
       * @param max_count maximum path count
       */
      void findAllPaths(Size source, Size sink, std::vector<std::vector<Size>>& all_paths, Size max_count)
      {
        std::unordered_set<Size> visited;
        std::vector<Size> path;

        // the adjacency is built once (at the first search that finds it empty): every edge once, by source and then target, so the
        // targets of a vertex are in ascending order
        if (! adjacency_used_)
        {
          sources_.clear();
          offsets_.clear();
          targets_.clear();
          auto add = [this](Size v1, Size v2) {
            if (sources_.empty() || sources_.back() != v1)
            {
              sources_.push_back(v1);
              offsets_.push_back(targets_.size());
            }
            targets_.push_back(v2);
          };
          // vertices are below vertex_count_ <= 2^bits: if an edge fits in 64 bits as (vertex1 << bits) | vertex2, the edges are sorted
          // as these numbers, which order them as the pairs do
          int bits = 1;
          while (bits < 64 && (Size(1) << bits) < vertex_count_) bits++;
          if (2 * bits <= 64)
          {
            std::vector<uint64_t> keys(edges_.size());
            for (Size e = 0; e < edges_.size(); e++) keys[e] = (uint64_t(edges_[e].first) << bits) | uint64_t(edges_[e].second);
            // not needed any more: an edge makes the adjacency non-empty, so it is never built again (see adjacency_used_)
            std::vector<std::pair<Size, Size>>().swap(edges_);
            sortKeys_(keys, 2 * bits);
            keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
            const uint64_t mask = (uint64_t(1) << bits) - 1;
            for (uint64_t key : keys) add((Size)(key >> bits), (Size)(key & mask));
          }
          else
          {
            std::sort(edges_.begin(), edges_.end());
            edges_.erase(std::unique(edges_.begin(), edges_.end()), edges_.end());
            for (const auto& [v1, v2] : edges_) add(v1, v2);
          }
          offsets_.push_back(targets_.size());
          adjacency_used_ = ! sources_.empty();
        }
        findAllPaths_(source, sink, visited, path, all_paths, max_count); // reverse traveling
      }

    private:
      Size vertex_count_;
      // 0, 1, 2, ... ,vertex_count - 1
      /// the edges as added (vertex1, vertex2), duplicates included
      std::vector<std::pair<Size, Size>> edges_;
      /// the adjacency of findAllPaths: the vertices with edges (ascending), and per vertex its targets (ascending) in targets_
      std::vector<Size> sources_, offsets_, targets_;
      /// whether the adjacency is built and kept (a search visited a vertex's adjacency, or it has an edge)
      bool adjacency_used_ = false;

      /// sort @p keys (all below 2^@p key_bits) ascending: large ones by an LSD radix sort over bytes (skipping bytes equal in all keys)
      static void sortKeys_(std::vector<uint64_t>& keys, int key_bits)
      {
        if (keys.size() < 16384)
        {
          std::sort(keys.begin(), keys.end());
          return;
        }
        const int passes = (key_bits + 7) / 8;
        std::vector<Size> counts(passes * 256, 0);
        for (uint64_t key : keys)
        {
          for (int p = 0; p < passes; p++) counts[p * 256 + ((key >> (8 * p)) & 255)]++;
        }
        std::vector<uint64_t> buffer(keys.size());
        for (int p = 0; p < passes; p++)
        {
          Size* count = counts.data() + p * 256;
          if (count[(keys[0] >> (8 * p)) & 255] == keys.size()) continue; // the same byte in all keys
          Size sum = 0;
          for (int d = 0; d < 256; d++)
          {
            const Size n = count[d];
            count[d] = sum;
            sum += n;
          }
          for (uint64_t key : keys) buffer[count[(key >> (8 * p)) & 255]++] = key;
          keys.swap(buffer);
        }
      }
      /// recursive part for the findAllPaths(..)
      void findAllPaths_(Size current,
                         Size destination,
                         std::unordered_set<Size>& visited,
                         std::vector<Size>& path,
                         std::vector<std::vector<Size>>& all_paths,
                         Size max_count)
      {
        if (max_count > 0 && all_paths.size() >= max_count) return;

        visited.insert(current);
        path.push_back(current);

        if (current == destination) { all_paths.push_back(path);}
        else
        {
          adjacency_used_ = true;
          const auto it = std::lower_bound(sources_.begin(), sources_.end(), current);
          if (it != sources_.end() && *it == current)
          {
            const Size k = (Size)(it - sources_.begin());
            for (Size e = offsets_[k]; e < offsets_[k + 1]; e++)
            {
              const Size i = targets_[e];
              if (visited.find(i) == visited.end()) { findAllPaths_(i, destination, visited, path, all_paths, max_count); }
            }
          }
        }

        // Backtrack
        path.pop_back();
        visited.erase(current);
      }
    };

    /**
       @brief calculate log mzs from mzs
       @param mz mz
       @param positive determines the charge carrier mass
     */
    static double getLogMz(double mz, bool positive);

    /**
       @brief get charge carrier mass : positive mode mass of (Constants\::PROTON_MASS_U) and negative mode mass of (-Constants\::PROTON_MASS_U)
       @param positive_ioniziation_mode Determines the charge carrier mass (true = positive or false = negative)
    */
    static float getChargeMass(bool positive_ioniziation_mode);

    /**
       @brief The modification a user names (FLASHTnT var_mod and fixed_mod): a name with its site as in Unimod/PSI-MOD, e.g. "Oxidation (M)",
       or a Unimod accession number with its site, e.g. "35 (M)", "1 (N-term)" ("UNIMOD:35 (M)" is the same; sites: a residue letter, N-term,
       C-term, Protein N-term, Protein C-term). A bare accession ("35") is taken only if the database knows one site for it.
       @throw Exception::ElementNotFound for an unknown name, or an accession unknown at the site
       @throw Exception::InvalidValue for a bare accession of several sites (the message lists them) or a malformed accession
    */
    static const ResidueModification* resolveModification(const std::string& spec);

    /**
       @brief A variable modification a user names (FLASHTnT var_mod): as resolveModification, and also a modification that the database has
       at any N-terminus only, restricted to the protein N-terminus by the site "Protein N-term <residue>", e.g. "Gln->pyro-Glu (Protein
       N-term Q)" for Gln->pyro-Glu (N-term Q); @p protein_n_term tells which.
       @throw as resolveModification for the name as given
    */
    static const ResidueModification* resolveVariableModification(const std::string& spec, bool& protein_n_term);

    /**
      @brief Read the spectra of an input run: an mzML file, a Thermo .raw file, or a Bruker timsTOF .d directory.

      An mzML file is read as MzMLFile reads it, with @p options. A .raw file (recognised by its extension) is read with Thermo's
      peak picking, so its spectra are centroided as ThermoRawFileParser and FileConverter convert them by default; the MS levels and
      the retention time range of @p options are applied while reading, the m/z and intensity ranges after it. Reading a .raw file
      needs OpenMS built with WITH_THERMO_RAW and the .NET 8 runtime; otherwise it stops with a message to convert the file to mzML.

      A .d directory or a zipped one (.d.zip; recognised by its name, see isBrukerTims) is read with OpenMS's timsTOF reader
      (BrukerTimsFile, DDA-PASEF or DIA-PASEF detected). Every spectrum is centroided over the ion mobility axis: a DDA MSn spectrum is
      one precursor summed over its frames and ion mobility scans and centroided in TOF space by the reader; an MS1 frame (and a DIA
      MSn frame) is summed over its ion mobility scans and centroided here in the same way. The ion mobility arrays are not kept. The
      spectra are numbered by their order in retention time (native IDs scan=1, 2, ...; the reader's native ID is kept as the meta value
      bruker_native_id), and @p options are applied after the numbering, so that the scan numbers are the same for every reading.
      Needs OpenMS built with WITH_OPENTIMS (the default); otherwise it stops with a message to convert the data to mzML.
    */
    static void loadExperiment(const std::string& file, PeakMap& map, const PeakFileOptions& options, ProgressLogger::LogType log);

    /// true if @p file is a Bruker timsTOF .d directory or a zipped one (.d.zip), by its name (a trailing separator allowed)
    static bool isBrukerTims(const std::string& file);

    /**
      @brief Fragmentation methods and the ion types the top-down tools search for them (the -fragmentation option of FLASHTnT and
      FLASHAnnotator, used when -ion_type is auto).

      fragmentationMethods: the method names the option takes after 'auto': HCD, CID, ETD, ECD, EThcD, EChcD, ECD-HCD, UVPD.
      ionTypesOf: the ion types of a method: HCD, CID: b y; EThcD, EChcD, ECD-HCD: b c y z zp1; ETD, ECD: c zp1 y; UVPD: a b c x y z;
      any other name (an unknown method): b y.
      fragmentationOf: the method of a spectrum from the activation of its first precursor as the file records it (mzML activation terms;
      a Thermo .raw file's activation type, its supplemental activation as EThcD/ETciD): ETD or ECD with a collisional activation (or
      EThcD, ETciD) is EThcD or ECD-HCD, ETD and ECD alone are ETD and ECD, beam-type CID is HCD, the other collisional activations
      (CID, trap-type, high- and low-energy CID, PQD) are CID, photodissociation (and a Thermo UVPD activation) is UVPD; "" if there is
      no precursor or the activation is none of these. Add-on ExD cells are recorded as the collisional activation they replace (HCD).
      countFragmentation: the MSn spectra of @p map per method ("" for unknown).
      resolveIonTypes: the ion types to search. @p ion_types other than {"auto"} are returned as given. With {"auto"}, a @p fragmentation
      other than "auto" gives its ion types; "auto" gives the ion types of the methods counted in @p counts (unknown as b y, with a
      warning naming the count); if they need different ion types, the union of them (with a warning), or with @p allow_mixed {}, and the
      caller resolves each spectrum. @p summary describes the choice for the log.
    */
    static const std::vector<std::string>& fragmentationMethods();
    static std::vector<std::string> ionTypesOf(const std::string& method);
    static std::string fragmentationOf(const MSSpectrum& spectrum);
    static std::map<std::string, Size> countFragmentation(const MSExperiment& map);
    static std::vector<std::string> resolveIonTypes(const std::vector<std::string>& ion_types, const std::string& fragmentation,
                                                    const std::map<std::string, Size>& counts, bool allow_mixed, std::string& summary);

    /**
      @brief A modification of a proteoform as the top-down tools read and write it in ProForma 2.0 (HUPO-PSI): its site, its residues and
      its mass. Parsing and writing are OpenMS::ProForma's (getProFormaModifications, toPeptidoform); the masses are computed here.

      What is read (getProFormaModifications) and written (toPeptidoform):
      - residues: the 20 standard amino acids, U and O
      - modification tags: a name looked up in Unimod, then PSI-MOD ([Oxidation], [L-methionine sulfoxide]); a name with its source
        ([U:Oxidation], [M:L-methionine sulfoxide]); an accession ([UNIMOD:35], [MOD:00719]); a formula ([Formula:C2H2O]); a glycan
        composition ([Glycan:Hex1HexNAc2]); a mass ([+15.9949], [Obs:+15.995]). Alternatives are separated by | and the first one that
        gives a mass is used ([Phospho|+79.966]); INFO: parts are ignored. A localisation group or score suffix (#g1, (0.9)) is dropped:
        the modification stays where its tag is written.
      - on a residue: PEM[Oxidation]K, several on one residue: S[Phospho][+1.0]; on a range: P(EMK)[+15.995]; on a terminus:
        [Acetyl]-PEPTIDE, PEPTIDE-[Amidated]; at an unknown position: [Phospho]?PEPTIDE, [Phospho]^2?PEPTIDE; labile:
        {Glycan:Hex}PEPTIDE; global fixed: <[Carbamidomethyl]@C>PEPTCDE; charge: PEPTIDE/2 (read only)
      - written in the order of ProForma 2.0: global, unknown position, labile, N-terminal, sequence, C-terminal

      Not supported (Exception::ParseError): chimeras (+), several chains (//), cross-links and branches (#XL, #BRANCH), isotope labels
      (<13C>), adducts as charge (/[Na:z+1]), names ((>name)), ambiguous residues ((?DQ)) and residues other than the 20, U and O.
    */
    struct OPENMS_DLLAPI ProFormaModification
    {
      /// where a modification sits
      enum class Site
      {
        RESIDUE, ///< on residues start..end (a range if end > start)
        N_TERM,
        C_TERM,
        UNKNOWN, ///< at an unknown position of the sequence
        LABILE   ///< lost on fragmentation: in the precursor mass only
      };

      Site site = Site::RESIDUE;
      int start = 0;       ///< first residue, 0-based (Site::RESIDUE)
      int end = 0;         ///< last residue, 0-based (Site::RESIDUE)
      double mass = 0;     ///< monoisotopic mass shift in Da
      std::string tag;     ///< the tag content as written (name, accession, formula, ...); empty: the mass is written (formatProFormaMass)
      bool global = false; ///< from a global fixed modification <[tag]@residues> (written there, not on the residue)
    };

    /**
      @brief The sequence and the modifications of a ProForma string read by OpenMS::ProForma::parseIon (see ProFormaModification).

      All modifications with their masses, in the order global (per residue), unknown position (a count ^n repeats them), labile,
      N-terminal, residues and ranges (as written; a range after the residues inside it), C-terminal. A terminal modification has the site
      N_TERM or C_TERM; a residue tag without mass ([#g1], [INFO:x]) is left out. The charge, if any, is in @p ion.
      @throw Exception::ParseError for a feature not supported (see ProFormaModification) or a modification without mass
    */
    static std::vector<ProFormaModification> getProFormaModifications(const ProForma::PeptidoformIon& ion, std::string& sequence);

    /**
      @brief A ProForma string read by OpenMS::ProForma::parseIon, after dropping the alternative "-" of a tag ("+42.010565|-" becomes
      "+42.010565"): FLASHTnT writes a fixed modification without an accession that way, and ProForma has no "-" tag.
      @throw Exception::ParseError as OpenMS::ProForma::parseIon
    */
    static ProForma::PeptidoformIon parseProForma(const std::string& text);

    /**
      @brief The proteoform @p sequence with @p modifications as a ProForma::Peptidoform, to be written by ProForma::toString (lossless):
      each tag as written (a mass formatted by formatProFormaMass), modifications of one residue or range in their order here, global ones
      as one <[tag]@residues> per tag.
      @throw Exception::ParseError for a residue other than the 20 standard ones, U and O
      @throw Exception::InvalidValue for residues outside the sequence, or a range that overlaps another (ProForma ranges are the same or
      disjoint); several modifications on one residue are allowed
    */
    static ProForma::Peptidoform toPeptidoform(const std::string& sequence, const std::vector<ProFormaModification>& modifications);

    /**
      @brief The mass of a ProForma tag content (see ProFormaModification), read by OpenMS::ProForma, e.g. Oxidation, U:Oxidation,
      UNIMOD:35, MOD:00719, Formula:O, Glycan:Hex1HexNAc1, INFO:x|Phospho|+79.97, Obs:+79.97, Cation:Ca[II]; @p residue (one-letter code or
      0) narrows a name or accession lookup. A name or accession that OpenMS::ProForma does not read (a Unimod name starting with a digit,
      2-methylbutyryl; OpenMS's UniMod:35) is looked up in ModificationsDB as it is.
      @throw Exception::ParseError if no alternative gives a mass
    */
    static double getProFormaTagMass(const std::string& tag, char residue = 0);

    /// a mass as written in a ProForma tag: sign, at most six decimals, no trailing zeros (+15.994915, -18.010565, +1.0); six decimals keep
    /// a written and read-back PrSM identical for annotation (four moved the mass calibration of FLASHAnnotator in 32 of 4245 PrSMs)
    static std::string formatProFormaMass(double mass);
  };

  // Hot in every peak-group sort and mass computation: defined inline (same arithmetic as before).
  inline float FLASHHelperClasses::getChargeMass(const bool positive_ioniziation_mode)
  {
    return (float)(positive_ioniziation_mode ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U);
  }

  inline double FLASHHelperClasses::LogMzPeak::getUnchargedMass() const
  {
    if (abs_charge == 0)
    {
      return .0;
    }
    if (mass <= 0)
    {
      return (mz - getChargeMass(is_positive)) * (float)abs_charge;
    }
    return mass;
  }

  inline bool FLASHHelperClasses::LogMzPeak::operator<(const LogMzPeak& a) const
  {
    if (this->logMz == a.logMz)
    {
      return this->intensity < a.intensity;
    }
    return this->logMz < a.logMz;
  }

  inline bool FLASHHelperClasses::LogMzPeak::operator>(const LogMzPeak& a) const
  {
    if (this->logMz == a.logMz)
    {
      return this->intensity > a.intensity;
    }
    return this->logMz > a.logMz;
  }

  inline bool FLASHHelperClasses::LogMzPeak::operator==(const LogMzPeak& a) const
  {
    return this->logMz == a.logMz && this->intensity == a.intensity;
  }
} // namespace OpenMS