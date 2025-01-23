// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <tuple>

namespace OpenMS
{
  namespace Internal
  {
    // The general class template
    template <size_t param_index, size_t grid_size, typename EvalResult, typename Tuple, typename... TupleTypes>
    struct Looper
    {
    };

    // Specialization for the base case
    // - shape_index == shape_size
    // - TupleTypes is empty here
    // - All indices in Functor are bound (i.e. can be called with empty arguments)
    template <size_t grid_size, typename EvalResult, typename Tuple, typename... TupleTypes>
    struct Looper<grid_size, grid_size, EvalResult, Tuple, TupleTypes...>
    {
      template <typename Functor>
      double operator()(const Tuple&, Functor functor, EvalResult /*bestValue*/, std::array<size_t, grid_size>& /*bestIndices*/)
      {
        return functor();
      }
    };

    // Specialization for the loop case
    // - increment shape_index
    // - create new Functor with one argument less and the first being bound to it
    // - loop over all values in the current vector and update best score and best indices
    template <size_t param_index, size_t grid_size, typename EvalResult, typename Tuple, typename FirstTupleType, typename... TupleTypes>
    struct Looper<param_index, grid_size, EvalResult, Tuple, FirstTupleType, TupleTypes...>
    {
      template <typename Functor>
      EvalResult operator()(const Tuple& grid, Functor functor, EvalResult bestValue, std::array<size_t, grid_size>& bestIndices)
      {
        for (size_t index = 0; index < std::get<param_index>(grid).size(); ++index)
        {
          double currVal = Looper<param_index + 1, grid_size, EvalResult, Tuple, TupleTypes...>()
              (
                  grid,
                  [&grid, index, &functor](TupleTypes... rest){ return functor(std::get<param_index>(grid)[index], rest...);},
                  bestValue,
                  bestIndices
              );

          if ( currVal > bestValue )
          {
            bestValue = currVal;
            bestIndices[param_index] = index;
          }
        }
        return bestValue;
      }
    };
  } // namespace Internal

  template <typename... TupleTypes>
  class GridSearch
  {
  public:
    explicit GridSearch(std::vector<TupleTypes>... gridValues):
        grid_(std::make_tuple<std::vector<TupleTypes>...>(std::move(gridValues)...))
    {}

    //Specific implementation for function objects
    template <typename Functor>
    typename std::result_of<Functor(TupleTypes...)>::type evaluate(Functor evaluator,
                                                                   typename std::result_of<Functor(TupleTypes...)>::type startValue,
                                                                   std::array<size_t,std::tuple_size<std::tuple<std::vector<TupleTypes>...>>::value>& resultIndices)
    {
      return Internal::Looper<0,
          std::tuple_size<std::tuple<std::vector<TupleTypes>...>>::value,
          typename std::result_of<Functor(TupleTypes...)>::type,
          std::tuple<std::vector<TupleTypes>...>,
          TupleTypes...> ()
          (grid_, evaluator, startValue, resultIndices);
    }


    //Specific implementation for function pointers
    template <typename EvalResult>
    EvalResult evaluate(EvalResult evaluator(TupleTypes...),
                        EvalResult startValue,
                        std::array<size_t,std::tuple_size<std::tuple<std::vector<TupleTypes>...>>::value>& resultIndices)
    {
      return Internal::Looper<0,
          std::tuple_size<std::tuple<std::vector<TupleTypes>...>>::value,
          EvalResult,
          std::tuple<std::vector<TupleTypes>...>,
          TupleTypes...>()
          (grid_, evaluator, startValue, resultIndices);
    }


    unsigned int getNrCombos()
    {
      if (combos_ready_)
      {
        return combos_;
      }
      else
      {
        return nrCombos();
      }
    }

  private:
    std::tuple<std::vector<TupleTypes>...> grid_;
    unsigned int combos_ = 1;
    bool combos_ready_ = false;

    template<std::size_t I = 0>
    typename std::enable_if<I == sizeof...(TupleTypes), unsigned int>::type
    nrCombos()
    {
      combos_ready_ = true;
      return combos_;
    }

    template<std::size_t I = 0>
    typename std::enable_if<I < sizeof...(TupleTypes), unsigned int>::type
    nrCombos()
    {
      combos_ *= std::get<I>(grid_).size();
      return nrCombos<I + 1>();
    }
  };
} // namespace OpenMS

