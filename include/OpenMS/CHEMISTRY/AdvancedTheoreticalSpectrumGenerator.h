// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_CHEMISTRY_ADVANCEDTHEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_ADVANCEDTHEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/SIMULATION/SimTypes.h>


namespace OpenMS
{
  /**
   @brief Generates theoretical spectra according to a probabilistic model. 

   The models are generated with the @ref UTILS_SpectrumGeneratorNetworkTrainer

   @htmlinclude OpenMS_AdvancedTheoreticalSpectrumGenerator.parameters

   @ingroup Chemistry
   */
  class OPENMS_DLLAPI AdvancedTheoreticalSpectrumGenerator: public TheoreticalSpectrumGenerator
  {
    public:
      ///Function object to convert the indices of the internal arrays
      class IndexConverter
      {
        public:
          UInt operator()(const UInt &type_id_a, const UInt &intensity_level_a, const UInt &intensity_level_parent, const UInt &number_intensity_levels)
          {
            return (type_id_a * number_intensity_levels * number_intensity_levels +
                intensity_level_parent * number_intensity_levels+
                intensity_level_a);
          }
      };

    public:
      /**
       @brief Nested class represents the TAN used for the probabilistic network
       */
      class OPENMS_DLLAPI TreeAugmentedNetwork
      {
        public:

          /**
           @name Small nested container representing an edge in a TreeAugmentedNetwork
           */
          //@{
          ///Representation of an edge in the probabilistic network with two end nodes and a weight
          struct TanEdge
          {
            UInt left_node;
            UInt right_node;
            DoubleReal score;

            //Overloaded operator < for the TanEdge
            bool operator<(const TanEdge &rhs) const
            {
              return score < rhs.score;
            }
          };
          //@}

          /** @name Constructors and Destructors*/
          //@{
          TreeAugmentedNetwork() :
            edges_(0),
            nodes_in_dfs_order_(0)
          {
          }


          TreeAugmentedNetwork(std::vector<TanEdge> edges_in) :
            edges_(edges_in),
            nodes_in_dfs_order_(0)
          {
          }

          ///Destructor
          virtual ~TreeAugmentedNetwork()
          {
          }

          ///Copy constructor
          TreeAugmentedNetwork(const TreeAugmentedNetwork &t_aug) :
            edges_(t_aug.edges_),
            nodes_in_dfs_order_(t_aug.nodes_in_dfs_order_)
          {
          }

          ///Assigment operator
          TreeAugmentedNetwork & operator = (const TreeAugmentedNetwork &t_aug)
          {
            if(this != &t_aug)
            {
              edges_ = t_aug.edges_;
              nodes_in_dfs_order_ = t_aug.nodes_in_dfs_order_;
            }
            return *this;
          }

          //@}

          /**
           @brief Computes the minimal spanning tree in a Kruskal like algorithm.

           Since the input Graph is very small a very simplified version of the algorithm is used
           @param tree_structure contains the resulting tree with a is-child-of entry for each node (-1 for root node)
         */
          void generateTree(std::vector<Int> &tree_structure);

          /**
           * @brief Return the node indices ordered according to their discovery in a DepthFirstSearch
           *
           * In this order the nodes can then be used in the simulation process, since then every type will be generated after the
           * one it conditionally depends on
           */
          void getDFSOrder(std::vector<UInt> &ordered_nodes) const
          {
            ordered_nodes=nodes_in_dfs_order_;
          }

          /** @name member variables
           */
          //@{
          ///Vecor of edges
          std::vector<TanEdge> edges_;
          ///Nodes in DFS order
          std::vector<UInt>nodes_in_dfs_order_;
          //@}
      };//End of TreeAugmentedNetwork


      /** @name Constructors and Destructors
       */
      //@{
      /// Default constructor
      AdvancedTheoreticalSpectrumGenerator();

      /// Copy constructor
      AdvancedTheoreticalSpectrumGenerator(const AdvancedTheoreticalSpectrumGenerator& source);

      /// Destructor
      virtual ~AdvancedTheoreticalSpectrumGenerator();
      //@}

      /// Assignment operator
      AdvancedTheoreticalSpectrumGenerator& operator =(const AdvancedTheoreticalSpectrumGenerator& tsg);

      /// Generate the MS/MS according to the given probabilistic model
      void simulate(RichPeakSpectrum &spectrum, const AASequence &peptide, const gsl_rng *rng, Int charge = 1);

      ///Load the probabilistic model from file
      void loadProbabilisticModel();
/*
      //for test reasons
      void writeProbabilisticModel(const String &file_name);
*/

      /**
       @brief nested class

      */
      //@{
      ///IonType is defined by a ResidueType, a neutral loss and charge state
      struct IonType{
        Residue::ResidueType residue;
        EmpiricalFormula loss;
        Int charge;

        /** @name Constructors and Destructors*/
        //@{
        //Default constructor
        IonType():
          residue((Residue::ResidueType) 0),
          loss(),
          charge(0)
          {
          }

        //Custom construtor
        IonType(Residue::ResidueType residue, EmpiricalFormula loss, Int charge):
          residue(residue),
          loss(loss),
          charge(charge)
          {
          }

        //Copy constructor
        IonType(const IonType & rhs):
            residue(rhs.residue),
            loss(rhs.loss),
            charge(rhs.charge)
            {
            }

        //Assignment operator
        IonType & operator = (const IonType & rhs)
        {
          if(this != &rhs)
          {
            residue = rhs.residue;
            loss = rhs.loss;
            charge = rhs.charge;
          }
          return *this;
        }
        //@}

        bool operator < (const IonType &rhs)
        {
          if(residue != rhs.residue)
            return residue<rhs.residue;
          else if(loss.getString()!=rhs.loss.getString())
            return loss.getString()<rhs.loss.getString();
          else
            return charge <rhs.charge;
        }
      };
      //@}

    private:
      ///Returns the ResidueType (e.g. AIon, BIon) as string for peak annotation
      String ResidueTypeToString_(Residue::ResidueType type);

      ///Vector of conditional probabilities for each sector
      std::vector<std::vector<DoubleReal> >conditional_probabilities_;

      ///The network models for each sector
      std::vector<TreeAugmentedNetwork> tan_;

      ///Number of discretized intensity levels;
      UInt number_of_intensity_levels_;

      ///Number of sectors for each Spectrum;
      UInt number_of_sectors_;

      ///The selected IonTypes
      std::vector<IonType>ion_types_;
  };

}
#endif
