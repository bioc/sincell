#include <Rcpp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pseudoreplicatesbynoise_cv2(NumericMatrix originaldata, int rows, int colums, NumericVector deciles, int lengthdeciles, NumericVector coorsorted, NumericVector vargenessorted, NumericVector means, int positive, int seed) {

  int x,y,i,range,start;
  float tempvar;
  NumericMatrix out(rows, colums);
  
  unsigned int iseed = (unsigned int)time(NULL);
  srand ((unsigned int)(iseed+seed));

    for(x=0; x<rows; x++){
      
      for(i=0; i<lengthdeciles-1; i++){
        if(coorsorted[x]<=deciles[i+1]){
          break;
        }
      }
      
      for(y=0; y<colums; y++){
      
        range = static_cast <int>(deciles[i+1]-deciles[i]);

        start = rand() % range + static_cast <int>(deciles[i]-1); 
        if(means[x]<0){tempvar=-means[x];} else {tempvar=means[x];}
        //tempvar = static_cast <float> (static_cast <float> (sqrt(vargenessorted[start]))*tempvar);
		tempvar = static_cast <float> (static_cast <float> (vargenessorted[start]*tempvar*tempvar));
  
        tempvar = static_cast <float> (3.0 * tempvar);
        tempvar = static_cast <float> (sqrt(tempvar));
                
        out(x,y) = originaldata(x,y) - tempvar + 2.0 * tempvar * static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));

        if(positive == 1){
          if(out(x,y)<0){
            out(x,y) = 0;
          }
        }
      
    }
    
  }

  return out;

}
