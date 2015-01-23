#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sstalgorithm(NumericVector membership, int num_cells, NumericMatrix distance){
  
  int x,y;
  int distx, disty;
  double distmin = 0;
  NumericVector out(3);
  
  for(x=0; x < (num_cells-1); x++){
    for(y=(x+1); y < num_cells; y++){
      if(membership[y]!=membership[x]){
              
        if(distance(x,y)< distmin || distmin == 0){
          distmin = distance(x,y);
          distx = x+1;
          disty = y+1;
        }
              
      }
    }
  }
  
  out[0] = distmin;
  out[1] = distx;
  out[2] = disty;
  
  return out;
  
}
