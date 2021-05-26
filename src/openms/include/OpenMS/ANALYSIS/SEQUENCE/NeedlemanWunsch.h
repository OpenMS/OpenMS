#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  class OPENMS_DLLAPI NeedlemanWunsch
  {

  public:
    enum class ScoringMatrix
    {
      identity,
      PAM30MS,
      SIZE_OF_SCORINGMATRIX
    };

    NeedlemanWunsch(ScoringMatrix matrix, int penalty);
    NeedlemanWunsch();

    ~NeedlemanWunsch()=default;

    static const std::vector<std::string> NamesOfScoringMatrices;

    int align(const String& seq1, const String& seq2);

    void setMatrix(const ScoringMatrix& matrix);
    void setMatrix(const std::string& matrix);

    void setPenalty(const int penalty);

    ScoringMatrix getMatrix() const;

    int getPenalty() const;

  private:
    unsigned seq1_len_ = 0;
    unsigned seq2_len_ = 0;
    int gap_penalty_ = 0;
    int my_matrix_ = 0;
    std::vector<int> first_row_{};
    std::vector<int> second_row_{};
  };
}