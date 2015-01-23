#include <Rcpp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix pseudoreplicatesbymodel(int rows, int colums, NumericVector alpha, NumericVector vargenes, NumericVector meangenes, int positive, Function f, int seed) {

  int x,y;
  NumericMatrix out(rows, colums);
  
    unsigned int iseed = (unsigned int)time(NULL);
    srand ((unsigned int)(iseed+seed));
    for(x=0; x<rows; x++){
  	for(y=0; y<colums; y++){
			if((static_cast <double> (alpha[x]) > (1/ static_cast <double> (colums)  ))){
				if(static_cast <double> (alpha[x]) > static_cast <float> (rand()) / static_cast <float> (RAND_MAX)){
						out(x,y) = as<double>(f(1,meangenes[x],sqrt(vargenes[x])));
					
				}else {
					out(x,y) = 0;
				}
			} else {
				out(x,y) = 0;
			}
			if(positive == 1){
				if(out(x,y)<0){
					out(x,y) = 0;
				}
			}
		}
	}
	return out;
}


