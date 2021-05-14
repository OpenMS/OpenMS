#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>


namespace OpenMS
{
  class OPENMS_DLLAPI NeedlemanWunsch
  {

  public:
    enum class ScoringMatrix
    {
      PAM30MSMatrix,
      identityMatrix
    };

    NeedlemanWunsch(ScoringMatrix matrix, int penalty);

    ~NeedlemanWunsch()=default;

    double align_(const String& seq1, const String& seq2);

    void setMatrix_(const ScoringMatrix& matrix);

    void setPenalty_(const int& penalty);

    ScoringMatrix getMatrix_() const;

    int getPenalty_() const;

  private:
    int getIndex_(const char& a, const char& b) const;
    unsigned seq1len_ = 0;
    unsigned seq2len_ = 0;
    int gapPenalty_ = 0;
    std::vector<int>* matrixPtr_ = nullptr;
  };
}