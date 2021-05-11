#include <OpenMS/ANALYSIS/SEQUENCE/NeedlemanWunsch.h>

using namespace std;
namespace OpenMS
{
  NeedlemanWunsch(ScoringMatrix matrix, int penalty)
{
  gapPenalty_ = penalty;
  switch(matrix)
  {
  case identity:
  matrixPtr_ = &adaptedIdentity;
  break;
  case PAM30MS:
  matrixPtr_ = &PAM30MS;
  break;
  default:
  String msg = "Matrix '" + matrix + "' is not known! Valid choices are: "
                                     "'identity', 'PAM30MS'.";
  throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
  msg);
  }
}

int NeedlemanWunsch::getIndex_(char& a const, char& b const) const //noch optimieren
{
vector<pair<char,int>> vec =
    {
        {'A', 0}, {'R', 1}, {'N', 2},
        {'D', 3}, {'C', 4}, {'Q', 5},
        {'E', 6}, {'G', 7}, {'H', 8},
        {'I', 9}, {'L', 10}, {'K', 11},
        {'M', 12}, {'F', 13}, {'P', 14},
        {'S', 15}, {'T', 16}, {'W', 17},
        {'Y', 18}, {'V', 19}, {'B', 20},
        {'Z', 21}, {'X', 22}, {'*', 23}
    };
int x = -1;
int y = -1;
for (int i = 0; i < vec.size(); ++i)
{
if (vec[i].first == a)
x = vec[i].second;
if (vec[i].first == b)
y = vec[i].second;
}
if (x == -1)
x = 23;
if (y == -1)
y = 23;
return x + y*vec.size();
}

int NeedlemanWunsch::align_(const String& seq1, const String& seq2)
{
  seq1len_ = seq1.length();
  seq2len_ = seq2.length();

  std::vector<int> matrix((seq1len_+1)*(seq2len_+1), 0)//matrix mit 0en initialisieren
  for (unsigned i = 1; i <= seq1len_; ++i) //vertikale mit gapkkosten initialisieren
    matrix[i*(seq2len_+1)]=i*gapPenalty_;
  for (unsigned i =0; i<=seq2len_;++i)//horizontale mit gapkosten initialieren
    matrix[i]=i*gapPenalty_;
  for (unsigned i=1;i<=seq1len_;++1)
  {
    for (unsigned j=1;j<=seq2len_;++j)
    {
      matrix[i*(seq2len_ +1)+j]=max(max((matrix[i*(seq2len_+1)+j-1]+gapPenalty_), (matrix[(i-1)*(seq2len_+1)+j]+gapPenalty)),
                                    (matrix[(i-1)*(seq2len_+1)+j-1])+ *matrixPtr_[getIndex_(seq1[i], seq2[j])])
    }
  }
  return matrix[(seq1len_+1)*(seq2len_+1)-1];
}
}