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
      PAM30MS,
      identity
    };

    NeedlemanWunsch(ScoringMatrix matrix, int penalty);

    ~NeedlemanWunsch()=default;

    int align(const String& seq1, const String& seq2);

    void setMatrix(const ScoringMatrix& matrix);

    void setPenalty(const int& penalty);

    ScoringMatrix getMatrix() const;

    int getPenalty() const;


  private:
    unsigned seq1len_ = 0;
    unsigned seq2len_ = 0;
    int gapPenalty_ = 0;
    int(* matrixPtr_)[26][26] = nullptr;
  };
}