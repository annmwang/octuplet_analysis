// pacman clustering
#ifndef pacman_h
#define pacman_h

#define NUMBOARD 8

using namespace std;

class pacman {
public:
  	void forward(vector< vector<double> >  mm_array[NUMBOARD], int mmhits[NUMBOARD], double clust_range, int tdo_cut, int iboard, double pad_at);
  	double backward(vector< vector<double> >  mm_array[NUMBOARD], int mmhits[NUMBOARD], double clust_range, int tdo_cut, int iboard, double pad_at);
 };

#endif
  // }

  // ~pacman(){}
  
  // returns charge in fC

	 //KEY//
	 //where a is board index, b loops through the hits
	 //mm_array[a][b][0] is CH #
	 //mm_array[a][b][1] is PDO
	 //mm_array[a][b][2] is TDO
	 //mm_array[a][b][3] is cluster number
inline void pacman::forward(vector< vector<double> >  mm_array[NUMBOARD], int mmhits[NUMBOARD], double clust_range, int tdo_cut, int iboard, double pad_at){
   	int nclus = 0;
   	int clust_start, clust_beg, clust_mult;
   	int clust_end = -1;
   	double clust_cha,clust_time;
   	int irep=0;
   	int k = iboard;
	clust_end = -1;
     	cout << "BOARD " << k << ", Numhits " << mmhits[k] << endl; 
     	for (int i=0; i<mmhits[k]; i++) {
       		if (i >clust_end && mm_array[k][i][1]>= 2.5  && (mm_array[k][i][2]-0.5*pad_at+120.) <= tdo_cut) { //WAS 5,then 10
		 		nclus++; // nclus indexes the found clusters
	 			clust_start=i;
	 			clust_end=i;
	 			clust_mult=1;
	 			mm_array[k][i][3] = nclus;
	 			clust_cha = mm_array[k][i][1];
	 			clust_time = mm_array[k][i][2];

	 			for(int j=i+1; j<mmhits[k]; j++) {
	   			//check for duplicate channels
	   				if (mm_array[k][j][0 ]== mm_array[k][i][0] && irep==0) {
	     				printf(" duplicate ch ev");
	     				//printf(" duplicate ch ev: %d\n ",mm_nev);
	     				irep=1;
	     
	     				for(int l=0; l<mmhits[k];l++) {
	       //UNCOMMENT THIS OUT FOR DUPLICATE ANALYSIS
	       // printf(" %5.1f %5.1f %5.1f  \n",
	       // 	      mm_array[k][l][0], mm_array[k][l][1], mm_array[k][l][2]);
	     				}		 
	   				}
	   
	   				if ( mm_array[k][j][1]>= 2.5 && (mm_array[k][j][2]- 0.5*pad_at + 120.) <= tdo_cut && 
						( mm_array[k][j][0]<=(mm_array[k][clust_end][0]+clust_range) )) {
	       				mm_array[k][j][3] = nclus; 
	       				clust_mult++;
	       				clust_cha = clust_cha+mm_array[k][j][1];
	       				clust_time = clust_time+mm_array[k][j][2];    
	       				clust_end = j;
	     			}	
	 			}
	 			mm_array[k][i][4] = clust_mult;
	 			mm_array[k][i][5] = clust_end;
	 			mm_array[k][i][6] = clust_cha;
	 			mm_array[k][i][7] = clust_time/float(clust_mult); 
       		}
     	
   	}
}
inline double pacman::backward(vector< vector<double> >  mm_array[NUMBOARD], int mmhits[NUMBOARD], double clust_range, int tdo_cut, int iboard, double pad_at){
  int nclus = 0;
  int clust_start, clust_beg, clust_mult;
  int clust_end = -1;
  int irep=0;
  int k = iboard;
  double run_clus_temp=0;
  double ind;
  int ind1_temp;
    for(int i=0; i<mmhits[k]; i++) {
      if(mm_array[k][i][3]> run_clus_temp) {
	clust_mult = mm_array[k][i][4];
	run_clus_temp = mm_array[k][i][3];
	ind = mm_array[k][i][0];
	ind1_temp = i;
	clust_beg = i;
	for(int j=i-1; j>-1 ; j--) {
	  if( ( mm_array[k][j][0]>=mm_array[k][clust_beg][0]- clust_range)
	      &&  mm_array[k][j][3]==0 && mm_array[k][j][1] >=2.5 &&
	      (mm_array[k][j][2]-0.5*pad_at+120.) <= tdo_cut ) {
	    mm_array[k][j][3]=run_clus_temp;
	    mm_array[k][j][4]= mm_array[k][ind1_temp][4]+1.;
	    mm_array[k][ind1_temp][4]=0.;
	    mm_array[k][j][5]= mm_array[k][ind1_temp][5];
	    mm_array[k][ind1_temp][5]=0;
	    mm_array[k][j][6]= mm_array[k][ind1_temp][6]+ mm_array[k][j][1];
	    mm_array[k][ind1_temp][6]=0.;
	    mm_array[k][j][7]= (mm_array[k][ind1_temp][7]* (mm_array[k][j][4]-1.)+
				mm_array[k][j][2])/mm_array[k][j][4];
	    mm_array[k][i][7]=0.;
	    ind1_temp=j;
	    clust_beg=j;
	  }
	}
      }
    }
    return run_clus_temp;
}
