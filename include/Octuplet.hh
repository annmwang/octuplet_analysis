// header file
#ifndef Octuplet_h
#define Octuplet_h
using namespace std;

class Octuplet {
public :
  static bool  compare(const vector<double> &, const vector<double>&);
  static bool  compare1(const vector<double> &, const vector<double>&);
  static bool  compare2(const vector<double> &, const vector<double>&);
  int getBoardIndex(int boardid, int run_num); 
};

#endif

inline bool Octuplet::compare2(const vector<double> &x, const vector<double> &y)
{
  return x[2] > y[2];
}
inline bool Octuplet::compare1(const vector<double> &x, const vector<double> &y)
{
  return x[0] > y[0];
}

inline bool  Octuplet::compare(const vector<double> &x, const vector<double> &y)
{
  return x[0] < y[0];
}

inline int Octuplet::getBoardIndex(int boardid, int run_num) {
  if (run_num == 3505){
    if (boardid == 110) return 0;
    else if (boardid == 109) return 1;
    else if (boardid == 101) return 2;
    else if (boardid == 103) return 3;
    else if (boardid == 106) return 4;
    else if (boardid == 102) return 5;
    else if (boardid == 107) return 6;
    else if (boardid == 105) return 7;
    else return -1;
  }
  if (run_num == (3508)){
    if (boardid == 111) return 0;
    else if (boardid == 116) return 1;
    else if (boardid == 101) return 2;
    else if (boardid == 109) return 3;
    else if (boardid == 112) return 4;
    else if (boardid == 102) return 5;
    else if (boardid == 107) return 6;
    else if (boardid == 105) return 7;
    else return -1;
  }
  else return -1;
}
