#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(y);
//  DATA_VECTOR(keep); //dataset.keep=keep; 
//  DATA_VECTOR_INDICATOR(keep, y); //dataset.keep=keep; 
  DATA_IVECTOR(iy);        
  DATA_IVECTOR(ia);        
  DATA_IVECTOR(is);
  DATA_IVECTOR(it);
  DATA_INTEGER(T);
  DATA_SCALAR(S);
  DATA_VECTOR(wt);
  DATA_VECTOR_INDICATOR(keep,y);
  DATA_INTEGER(adr_yhat);             
  PARAMETER(beta); 
  PARAMETER_VECTOR(log_sd_err);
  PARAMETER(log_sd_N);  
  PARAMETER_MATRIX(q); 
  PARAMETER_VECTOR(Ny);
  // PARAMETER(c_adj);
   
  Type nll = Type(0.0);   
  
  vector<Type> sd_err = exp(log_sd_err);  
  Type sd_N = exp(log_sd_N);
  
  vector<Type> Np(T);
  Np(0)=0.0;
  for(int i = 1;i < T;++i){Np(i)=Ny(i-1);} 
     
  vector<Type> Nay = beta + Np(iy);
  int n = y.size();   
  vector<Type> mu(n);
  for(int i = 0;i < n;++i){
    if(it(i) == 0){
      mu(i) = Nay(i) + q(is(i),ia(i));
    }else{
      mu(i) = S*(Nay(i) + q(is(i),ia(i)));
    }
    nll -= keep(i)*wt(i)*dnorm(y(i), mu(i), sd_err(is(i)), true);
  }

  

  for(int i = 1;i < T;++i){
    nll -= dnorm(Np(i), Np(i-1), sd_N, true);
  }
  
  REPORT(beta); 
  REPORT(q);    
  REPORT(Ny);   
  REPORT(Np);    
  REPORT(sd_N);
  REPORT(sd_err); 
  REPORT(Nay); 
  REPORT(wt);   
  REPORT(mu);   
  ADREPORT(Ny); 
  if(adr_yhat==1){ADREPORT(mu);}    
  return nll;
}

