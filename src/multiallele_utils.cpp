#include <Rcpp.h>
using namespace Rcpp;

// Utilities for working with multiallelic genotypes.

// Compiled version of the dmultinom function from R, using some of the same
// source code.  Performs much less error checking than the R version.
// [[Rcpp::export]]
double dmultinom(NumericVector x, NumericVector prob){
  double out;
  double logout;
  LogicalVector nonzero = prob > 0;
  x = x[nonzero];
  prob = prob[nonzero];
  
  logout = lgamma(sum(x) + 1) + sum(x * log(prob) - lgamma(x + 1));
  out = exp(logout);
  
  return out;
}

// Probability distrubution under the Dirichlet multinomial.
// For estimating genotype likelihoods under overdispersion.
// [[Rcpp::export]]
double dDirichletMultinom(NumericVector x, NumericVector prob, double alpha){
  LogicalVector nonzero = prob > 0;
  x = x[nonzero];
  prob = prob[nonzero];
  double n = sum(x);
  double tot = lgamma(n + 1) + lgamma(alpha) - lgamma(n + alpha);
  NumericVector alphas = alpha * prob;
  double ind = sum(lgamma(x + 1) + lgamma(alphas) - lgamma(alphas + x));
  return exp(tot - ind);
}

// Function to get number of possible genotypes.
// (ploidy + nalleles - 1)!/(ploidy! * (nalleles - 1)!)
// [[Rcpp::export]]
int nGen(int ploidy, int nalleles) {
  // Adding 0.5 is necessary because of floating point math error and
  // truncation to integer.
  return exp(lgamma(ploidy + nalleles) - lgamma(ploidy + 1) - lgamma(nalleles)) + 0.5;
}

// Function to enumerate genotypes in the VCF order.
// Each row of output matrix is one genotype.
// The matrix has as many columns as the ploidy.
// Each cell contains the index of an allele that the genotype has.
// [[Rcpp::export]]
IntegerMatrix enumerateGenotypes(int ploidy, int nalleles){
  int ngen = nGen(ploidy, nalleles);
  IntegerMatrix out(ngen, ploidy);
  int startrow;
  IntegerMatrix submat;
  
  // a is allele, r is row, and c is column.
  for(int a = 0; a < nalleles; a++){
    if(ploidy == 1){
      out(a, 0) = a;
    } else {
      startrow = nGen(ploidy, a);
      submat = enumerateGenotypes(ploidy - 1, a + 1);
      for(int r = 0; r < submat.nrow(); r++){
        out(r + startrow, ploidy - 1) = a;
        for(int c = 0; c < ploidy - 1; c++){
          out(r + startrow, c) = submat(r, c);
        }
      }
    }
  }
  
  return out;
}

// Get the index of a given genotype, i.e. the row that the genotype would
// appear in, in the matrix output by enumerateGenotypes.
// [[Rcpp::export]]
int indexGenotype(IntegerVector genotype){
  int ploidy = genotype.size();
  int out = 0;
  int g;
  for(int m = 1; m < ploidy + 1; m++){
    g = genotype(m - 1);
    out += Rf_choose(g + m - 1, m);
  }
  return out;
}

// Get the copy number of each allele in a given genotype.
// Useful for converting genotypes exported by enumerateGenotypes to
// a format useful for multinomial probability calculations.
// [[Rcpp::export]]
IntegerVector alleleCopy(IntegerVector genotype, int nalleles){
  IntegerVector out (nalleles);
  int ploidy = genotype.size();
  int g;
  
  for(int i = 0; i < ploidy; i++){
    g = genotype(i);
    out(g)++;
  }
  
  return out;
}

// Internal, recursive function for generating gamete genotypes
IntegerMatrix makeGametesRecur(IntegerVector genotype, int gamploidy){
  int ploidy = genotype.size();
  int ngametes = Rf_choose(ploidy, gamploidy);
  IntegerMatrix out(ngametes, gamploidy);
  IntegerMatrix subout;
  int row = 0;
  
  if(gamploidy == 1){
    for(int i = 0; i < ngametes; i++){
      out(i, 0) = genotype(i);
    }
  } else {
    for(int i = 0; i < ploidy - gamploidy + 1; i++){
      subout = makeGametesRecur(genotype[Range(i + 1, (ploidy - 1))], gamploidy - 1);
      for(int r = 0; r < subout.nrow(); r++){
        out(row, 0) = genotype(i);
        for(int c = 0; c < subout.ncol(); c++){
          out(row, c + 1) = subout(r, c);
        }
        row++;
      }
    }
  }
  
  return out;
}

// For a given genotype, return a matrix containing all possible gamete
// genotypes.  Duplicates will be output reflecting relative gamete frequency.
// Ploidy should be even.
// [[Rcpp::export]]
IntegerMatrix makeGametes(IntegerVector genotype){
  int gamploidy = genotype.size() / 2;
  return makeGametesRecur(genotype, gamploidy);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
nGen(4, 3)
nGen(2, 2)
enumerateGenotypes(4, 3)
enumerateGenotypes(2, 2)
indexGenotype(c(0, 2))
indexGenotype(c(1, 1))
indexGenotype(c(0, 0, 2, 2))
alleleCopy(c(1, 1, 1,2), 4)

dmultinom(c(20, 25, 35), c(0.25, 0.25, 0.5))
dDirichletMultinom(c(20, 25, 35), c(0.25, 0.25, 0.5), 9)

makeGametes(c(0, 1))
makeGametes(0:3)
makeGametes(0:5)
*/
