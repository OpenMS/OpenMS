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
#include <concepts>
#include <ranges>
#include <functional>

namespace OpenMS
{
  namespace Internal
  {
    template<typename F, typename... Args>
    concept Evaluator = requires(F f, Args... args) {
      { std::invoke(f, args...) } -> std::convertible_to<double>;
    };

    // The general class template
    template <size_t param_index, size_t grid_size, typename EvalResult, typename Tuple, typename... TupleTypes>
    struct Looper
    {
    };

    // Specialization for the base case
    template <size_t grid_size, typename EvalResult, typename Tuple, typename... TupleTypes>
    struct Looper<grid_size, grid_size, EvalResult, Tuple, TupleTypes...>
    {
      template <typename Functor>
        requires Evaluator<Functor, TupleTypes...>
      constexpr auto operator()(const Tuple&, Functor functor, EvalResult /*bestValue*/, std::array<size_t, grid_size>& /*bestIndices*/) const
      {
        return std::invoke(functor);
      }
    };

    // Specialization for the loop case
    template <size_t param_index, size_t grid_size, typename EvalResult, typename Tuple, typename FirstTupleType, typename... TupleTypes>
    struct Looper<param_index, grid_size, EvalResult, Tuple, FirstTupleType, TupleTypes...>
    {
      template <typename Functor>
        requires Evaluator<Functor, FirstTupleType, TupleTypes...>
      constexpr auto operator()(const Tuple& grid, Functor functor, EvalResult bestValue, std::array<size_t, grid_size>& bestIndices) const
      {
        const auto& current_vector = std::get<param_index>(grid);
        
        for (size_t index = 0; index < current_vector.size(); ++index) {
          const auto& value = current_vector[index];
          auto currVal = Looper<param_index + 1, grid_size, EvalResult, Tuple, TupleTypes...>{}
              (
                  grid,
                  [&value, &functor](TupleTypes... rest){ 
                    return std::invoke(functor, value, rest...);
                  },
                  bestValue,
                  bestIndices
              );

          if (currVal > bestValue) {
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
    explicit GridSearch(std::vector<TupleTypes>... gridValues)
        : grid_(std::make_tuple<std::vector<TupleTypes>...>(std::move(gridValues)...))
    {}

    // Implementation for function objects using concepts
    template <typename Functor>
      requires Internal::Evaluator<Functor, TupleTypes...>
    constexpr auto evaluate(
        Functor evaluator,
        std::invoke_result_t<Functor, TupleTypes...> startValue,
        std::array<size_t, std::tuple_size_v<std::tuple<std::vector<TupleTypes>...>>>& resultIndices) const
    {
      return Internal::Looper<
          0,
          std::tuple_size_v<std::tuple<std::vector<TupleTypes>...>>,
          std::invoke_result_t<Functor, TupleTypes...>,
          std::tuple<std::vector<TupleTypes>...>,
          TupleTypes...>{}(grid_, evaluator, startValue, resultIndices);
    }

    // Implementation for function pointers using concepts
    template <typename EvalResult>
      requires std::convertible_to<EvalResult, double>
    [[nodiscard]] constexpr auto evaluate(
        EvalResult (*evaluator)(TupleTypes...),
        EvalResult startValue,
        std::array<size_t, std::tuple_size_v<std::tuple<std::vector<TupleTypes>...>>>& resultIndices) const
    {
      return Internal::Looper<
          0,
          std::tuple_size_v<std::tuple<std::vector<TupleTypes>...>>,
          EvalResult,
          std::tuple<std::vector<TupleTypes>...>,
          TupleTypes...>{}(grid_, evaluator, startValue, resultIndices);
    }

    [[nodiscard]] constexpr auto getNrCombos() const -> unsigned int
    {
      if (combos_ready_) {
        return combos_;
      }
      return calculateCombos();
    }

  private:
    std::tuple<std::vector<TupleTypes>...> grid_;
    mutable unsigned int combos_ = 1;
    mutable bool combos_ready_ = false;

    template<std::size_t I = 0>
    [[nodiscard]] constexpr unsigned int calculateCombos() const
    {
      if constexpr (I == sizeof...(TupleTypes)) {
        combos_ready_ = true;
        return combos_;
      } else {
        combos_ *= std::get<I>(grid_).size();
        return calculateCombos<I + 1>();
      }
    }
  };
} // namespace OpenMS
