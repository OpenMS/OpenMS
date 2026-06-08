#include "MonotoneRegressor.h"
#include "PAVARegressor.h"
#include "ISplineTRRRegressor.h"
#include <stdexcept>

namespace OpenMS { namespace Internal { namespace Percolator {

std::unique_ptr<MonotoneRegressor> make_monotone_regressor(RegressorType type) {
  std::unique_ptr<MonotoneRegressor> ptr;
  switch (type) {
    case RegressorType::PAVA:
      ptr = std::unique_ptr<MonotoneRegressor>(new PAVARegressor());
      break;
    case RegressorType::ISPLINE_TRR:
      ptr = std::unique_ptr<MonotoneRegressor>(new ISplineTRRRegressor());
      break;
    default:
      throw std::invalid_argument("make_monotone_regressor: unsupported RegressorType");
  }
  ptr->set_params(MonotoneParams());
  return ptr;
}

}}}  // namespace OpenMS::Internal::Percolator
