//scintillator analysis
#ifndef scint_h
#define scint_h

using namespace std;

class scint {
public :
  int pairs(vector< vector<double> > & hitarray, double low, double high);
};

#endif

//////////// look for W+E pairs
// pass an array, hitarray[a][b] where [a] indexes the counter hits,
// [a][0] is counter # (channel), [a][1] is TDC recorded
// low, high are the TDC cuts to get rid of crap events
// 
inline int scint::pairs(vector < vector<double> > & hitarray, double low, double high){
  int hits=0;
  int n = hitarray.size();
  for(int i=0; i<n;i++)   {
    for(int j=i+1; j<n;j++) {
      if(abs(hitarray[i][0]-hitarray[j][0])==6.) {
	if( (hitarray[i][1]>low && hitarray[i][1]<high)  &&
	    (hitarray[j][1]>low &&hitarray[j][1]<high) ) {
	  hitarray[i][2]=1.;
	  hitarray[j][2]=1.;
	  hits++;
	}
      }
    }
  }
  return hits;
}
