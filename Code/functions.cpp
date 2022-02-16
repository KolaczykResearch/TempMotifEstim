#include <Rcpp.h>
#include "graph_es.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <time.h>
#include <set>
#include <algorithm>
#include <map>
#include <random>
#include <sys/time.h>
#include <cmath>
using namespace std;
using namespace Rcpp;

// Exporting a C++ function to R. 
// [[Rcpp::export]]
NumericVector main2(NumericVector x) {
  Graph G;
  long a[3];
  G.Initialize();
  cout<<"Finish reading "<<G.getEn()<<" edges and the last edge is "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
  return x;
}

// [[Rcpp::export]]
double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

// [[Rcpp::export]]
NumericVector getVec(NumericVector y) {
  NumericVector xx(y.begin(), y.end());
  return xx;
}



// [[Rcpp::export]]
double get_wall_time(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// [[Rcpp::export]]
vector<double> estimate_motif_counts_var(NumericMatrix data, NumericMatrix motif, int delta, double p, int seed){
  // cout
  // cout<<"Running estimate_motif_counts() \n";
  // cout<<"delta = "<<delta<<endl;
  // cout<<"p = "<<p<<endl;
  
  //
  //==read the temporal graph data====================================================
  //
  Graph G;
  int n = data.nrow();
  for(int i=0; i<data.nrow(); i++){
    /*NumericVector a(data(i, _).begin(), data(i, _).end());
     if (a[0] != a[1]) {
     G.addEdge(a[0], a[1], a[2],0L);
     }*/
    if (data(i, 0) != data(i, 1)){
      G.addEdge(data(i, 0), data(i, 1), data(i, 2), 0L);
    } 
  }
  G.Initialize();
  /*for(int i=0;i<G.edges_.size();i++)
   {
   cout<<G.edges_[i].src<<" "<<G.edges_[i].dst<<" "<<G.edges_[i].tim<<" "<<G.edges_[i].id<<endl;
   }*/
  
  //uncomment
  // cout<<"Finish reading "<<G.getEn()<<" edges and the last edge is "<<data(n-1, 0)<<" "<<data(n-1, 1)<<" "<<data(n-1, 2)<<endl;
  
  // cout<<"Finish reading "<<G.getEn()<<" edges and the last edge is "<<a[0]<<" "<<a[1]<<" "<<a[2]<<endl;
  /*G.Printgraph();
   for (int i=0;i<G.getVn();i++)
   {
   sum_d=sum_d+G.nodeDegree(i);
   }
   cout<<"Sum of in and out degree: "<<sum_d<<endl;
   cout<<"the number of edges in hash table: "<<G.getHashn()<<endl;
   G.printGraph();*/
  
  //
  //==read a motif====================================================
  //
  G.cle_motif();
  for(int i=0; i<motif.nrow(); i++){
    if (motif(i, 0) != motif(i, 1)){
      // cout<<motif(i, 0)<<" "<<motif(i, 1)<<endl;
      G.motif(motif(i, 0), motif(i, 1));
    }else{cout<<"the motif has self-loops"<<motif(i, 0)<<" "<<motif(i, 1)<<endl; exit(1);}
  }
  // uncomment
  // cout<<"Finishing reading the motif"<<endl;
  
  //
  //==sampling====================================================
  //
  double rn;
  vector<long> ids;
  default_random_engine e(seed);
  uniform_real_distribution<double> u(0.0,1.0);
  for(long i=0;i<G.getEn();i++)
  {
    rn=u(e);
    //cout<<rn<<endl;
    if(rn<=p){
      /*if(G.edges_[i].tim==1197014240)
       {*/
      //cout<<i<<" ";
      ids.push_back(i);
    }
  }
  // uncomment
  // cout<<"Number of sampled edges:"<<ids.size()<<endl;
  
  //cout<<G.getHashn()<<endl;
  
  //
  //==estimate motif counts====================================================
  //
  double start,time0;
  G.setVm_();
  start=get_wall_time(); 
  //cout<<"begin"<<endl;
  // ExactCountMotifs() returns both (sum eta) and (sum eta^2)
  vector<double> output=G.ExactCountMotifs(delta,ids);
  double T=output[0]/p;
  double T2=output[1];
  time0=get_wall_time()-start;
  // cout<<"T="<<T<<" T/"<<G.edgesM_.size()<<"="<<T/G.edgesM_.size()<<" "<<time0<<endl;
  // out<<T/G.edgesM_.size()<<" "<<time0<<endl;
  double estimate = T/G.edgesM_.size(); 
  double sig_2_hat = T2*(1-p)/(pow(p,2)*pow(G.edgesM_.size(),2));
  // uncomment
  // cout<<"C_hat = "<< estimate<<endl;
  // cout<<"Time used = "<<time0<<endl;
  vector<double> return_vals;
  return_vals.push_back(estimate);
  return_vals.push_back(sig_2_hat);
  return return_vals;
  //return estimate;
}






