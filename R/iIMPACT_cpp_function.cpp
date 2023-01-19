#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
using namespace Rcpp;

arma::mat rdirichlet_cpp(NumericVector alpha, int size);
double lgamma_function(double alpha);
Rcpp::List iIMPACT(arma::mat V, arma::mat Y, arma::mat P,  int K, NumericVector e, double f, NumericVector eta, double tau, double alpha_gamma, double beta_gamma, NumericVector mu_initial, NumericVector sigma_initial, NumericVector alpha, NumericVector omega_initial, arma::mat Z_initial, NumericVector weight);
arma::mat rnorm_function(int n, NumericVector mu, arma::mat sigma);
double dnorm_function(arma::mat alpha, NumericVector mu, arma::mat sigma);
double lbeta_function(arma::mat alpha);
  


// [[Rcpp::export]]
Rcpp::List iIMPACT_Y(arma::mat Y, arma::mat P,  int K, NumericVector e, double f, NumericVector eta, double tau, double alpha_gamma, double beta_gamma, NumericVector mu_initial, NumericVector sigma_initial, arma::mat Z_initial) {
  // Read data information
  int N = P.n_rows;
  int n_neighbor = P.n_cols;
  int Y_dim = Y.n_cols;
  
  
  // Set algorithm settings
  int iter = 20000;
  int burn = iter/2;
  
  int i, q, qq, qqq, qqqq, l, h, m, count;
  int count_2 = 10;
  
  NumericMatrix mu_store(iter, Y_dim*K);
  NumericMatrix mu(K, Y_dim);
  NumericMatrix sigma(K, Y_dim);
  NumericMatrix sigma_store(iter, Y_dim*K);
  NumericMatrix Z_store(iter, N);
  NumericVector Z(N);
  
  NumericVector count_temp(K);
  arma::mat average_temp(K, Y_dim);
  
  NumericVector prob_temp(K);
  NumericVector neighbor_index(K);
  NumericVector sample_index(K);
  
  // Initialization
  for(q = 0; q < K; q++)
  {
    for (qq = 0; qq < Y_dim; qq++)
    {
      mu(q, qq) = mu_initial(qq);
      sigma(q, qq) = sigma_initial(qq);
    }
    sample_index(q) = q + 1; 
  }
  
  for (q = 0; q < N; q++)
  {
    Z(q) = Z_initial(q, 0);
  }
  
  
  // MCMC
  for(i = 0; i < iter; i++)
  {
    
    // reset to zero
    for(q = 0; q < K; q++)
    {
      count_temp(q) = 0;
      for (qq = 0; qq < Y_dim; qq++){
        average_temp(q, qq) = 0;
      }
    }
    
    // compute count and y_average
    for(q = 0; q < K; q++)
    {
      for (h = 0; h < N; h++){
        if (Z(h) == q + 1){
          count_temp(q) ++;
          for (qq = 0; qq < Y_dim; qq++){
          average_temp(q, qq) += Y(h, qq);
          }
        }
      }
    }
    
    for(q = 0; q < K; q++)
    {
    for (qq = 0; qq < Y_dim; qq++){
      if (count_temp(q) > 0){
      average_temp(q, qq) =  average_temp(q, qq)/count_temp(q);}
    }
    }
    
    // Update sigma
    for(q = 0; q < K; q++)
    {
      for (qq = 0; qq < Y_dim; qq++){
        double alpha_gamma_new = alpha_gamma + 0.5*count_temp(q);
        double sum_temp3 = 0;
        for (h = 0; h < N; h++){
          if (Z(h) == q + 1){
            sum_temp3 += (Y(h, qq) - average_temp(q, qq))*(Y(h, qq) - average_temp(q, qq)); 
          }
        }
        double beta_gamma_new = beta_gamma + 0.5*tau*count_temp(q)*(average_temp(q, qq) - eta(qq))*(average_temp(q, qq) - eta(qq))/(tau + count_temp(q)) + 0.5*sum_temp3;
       
       double tempp = rgamma(1, alpha_gamma_new, (1.0/beta_gamma_new))(0);
       sigma(q, qq) = 1.0/tempp;
      }
    }
    
    // Update mu
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        double eta_new = (eta(qq)*tau + count_temp(q)*average_temp(q, qq))/(count_temp(q) + tau);
        double sigma_for_mu = sigma(q, qq)/(tau + count_temp(q));
        mu(q, qq) = rnorm(1, eta_new, std::sqrt(sigma_for_mu))(0);
      }
    }
    
    // Update Z
    for (h = 0; h < N; h++){
      double sum_temp = 0;
      double sum_temp2 = 0;
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = 0;
        neighbor_index(q) = 0;
      }
      
      for (q = 0; q < n_neighbor; q++){
        if (P(h, q) != 0){
          neighbor_index(Z(P(h, q) - 1) - 1) = neighbor_index(Z(P(h, q) - 1) - 1) + 1;
        }
      }
      
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = e(q) + f*neighbor_index(q);
        for (qq = 0; qq < Y_dim; qq++){
          prob_temp(q) = prob_temp(q) -0.5*log(sigma(q, qq)) - 0.5*(Y(h, qq) - mu(q, qq))*(Y(h, qq) - mu(q, qq))/sigma(q, qq) ;
        }
        sum_temp = sum_temp + prob_temp(q); 
      }
      
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = prob_temp(q) - sum_temp/K;
        if (prob_temp(q) > 50){prob_temp(q) = 50;}
        if (prob_temp(q) < -50){prob_temp(q) = -50; }
        prob_temp(q) = exp(prob_temp(q));
      
        sum_temp2 = sum_temp2 + prob_temp(q);
      }
      
      if(sum_temp2 == 0){
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = 1.0/K;
        }
        Rcout <<i<< "% zero\n";
      }
      else{
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = prob_temp(q)/sum_temp2;
        }
      }
      
      Z(h) = sample(sample_index, 1, 1==0, prob_temp)(0); 
    }
    
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        mu_store(i, count) = mu(q, qq);
        sigma_store(i, count) = sigma(q, qq);
        count++;
      }
    }
    
    for(q = 0; q < N; q++)
    {
      Z_store(i, q) = Z(q);
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu_store, Rcpp::Named("sigma") = sigma_store, Rcpp::Named("Z") = Z_store);
}



// [[Rcpp::export]]
Rcpp::List iIMPACT(arma::mat V, arma::mat Y, arma::mat P,  int K, NumericVector e, NumericVector f, NumericVector eta, double tau, double alpha_gamma, double beta_gamma, NumericVector mu_initial, NumericVector sigma_initial, NumericVector alpha, NumericVector omega_initial, arma::mat Z_initial, NumericVector weight) {
  // Read data information
  int N = P.n_rows;
  int n_neighbor = P.n_cols;
  int Y_dim = Y.n_cols;
  int V_dim = V.n_cols;
  
  
  // Set algorithm settings
  int iter = 20000;
  int burn = iter/2;
  
  int i, q, qq, qqq, qqqq, l, h, m, count;
  int count_2 = 10;
  
  NumericMatrix mu_store(iter, Y_dim*K);
  NumericMatrix mu(K, Y_dim);
  NumericMatrix sigma(K, Y_dim);
  NumericMatrix sigma_store(iter, Y_dim*K);
  NumericMatrix Z_store(iter, N);
  NumericVector Z(N);
  NumericMatrix omega_store(iter, V_dim*K);
  NumericMatrix omega(K, V_dim);
  
  NumericVector count_temp(K);
  arma::mat average_temp(K, Y_dim);
  
  NumericVector V_temp(V_dim);
  arma::mat omega_temp(1, V_dim);
  NumericVector prob_temp(K);
  NumericVector neighbor_index(K);
  NumericVector sample_index(K);
  
  // Initialization
  for(q = 0; q < K; q++)
  {
    for (qq = 0; qq < Y_dim; qq++)
    {
      mu(q, qq) = mu_initial(qq);
      sigma(q, qq) = sigma_initial(qq);
    }
    sample_index(q) = q + 1; 
  }
  
  for (q = 0; q < N; q++)
  {
    Z(q) = Z_initial(q, 0);
  }
  
  for(q = 0; q < K; q++)
  {
    for (qq = 0; qq < V_dim; qq++)
    {
      omega(q, qq) = omega_initial(qq);
    }
    sample_index(q) = q + 1; 
  }
  
  
  
  // MCMC
  for(i = 0; i < iter; i++)
  {
    
    // reset to zero
    for(q = 0; q < K; q++)
    {
      count_temp(q) = 0;
      for (qq = 0; qq < Y_dim; qq++){
        average_temp(q, qq) = 0;
      }
    }
    
    // compute count and y_average
    for(q = 0; q < K; q++)
    {
      for (h = 0; h < N; h++){
        if (Z(h) == q + 1){
          count_temp(q) ++;
          for (qq = 0; qq < Y_dim; qq++){
            average_temp(q, qq) += Y(h, qq);
          }
        }
      }
    }
    
    for(q = 0; q < K; q++)
    {
      for (qq = 0; qq < Y_dim; qq++){
        if (count_temp(q) > 0){
        average_temp(q, qq) =  average_temp(q, qq)/count_temp(q);}
      }
    }
    
    // Update sigma
    for(q = 0; q < K; q++)
    {
      for (qq = 0; qq < Y_dim; qq++){
        double alpha_gamma_new = alpha_gamma + 0.5*count_temp(q);
        double sum_temp3 = 0;
        for (h = 0; h < N; h++){
          if (Z(h) == q + 1){
            sum_temp3 += (Y(h, qq) - average_temp(q, qq))*(Y(h, qq) - average_temp(q, qq)); 
          }
        }
        double beta_gamma_new = beta_gamma + 0.5*tau*count_temp(q)*(average_temp(q, qq) - eta(qq))*(average_temp(q, qq) - eta(qq))/(tau + count_temp(q)) + 0.5*sum_temp3;
        sigma(q, qq) = 1.0/rgamma(1, alpha_gamma_new, (1.0/beta_gamma_new))(0);
      }
    }
    
    // Update mu
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        double eta_new = (eta(qq)*tau + count_temp(q)*average_temp(q, qq))/(count_temp(q) + tau);
        double sigma_for_mu = sigma(q, qq)/(tau + count_temp(q));
        mu(q, qq) = rnorm(1, eta_new, std::sqrt(sigma_for_mu))(0);
      }
    }
    
    // Update omega
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < V_dim; qq++)
      {
        V_temp(qq) = alpha(qq);
      }
      
      for(qq = 0; qq < N; qq++)
      {
        if (Z(qq) == q + 1){
          for(qqq = 0; qqq < V_dim; qqq++){
            V_temp(qqq) = V_temp(qqq) + V(qq, qqq);
          }}
      }
      
      omega_temp = rdirichlet_cpp(V_temp, V_dim);
      
      
      for(qq = 0; qq < V_dim; qq++)
      {
        omega(q, qq) = omega_temp(0, qq);
        
      }
      
    }
    
    // Update Z
    for (h = 0; h < N; h++){
      double sum_temp = 0;
      double sum_temp2 = 0;
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = 0;
        neighbor_index(q) = 0;
      }
      
      for (q = 0; q < n_neighbor; q++){
        if (P(h, q) != 0){
          neighbor_index(Z(P(h, q) - 1) - 1) = neighbor_index(Z(P(h, q) - 1) - 1) + 1;
        }
      }
      
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = e(q) + f(h)*neighbor_index(q);
        for (qq = 0; qq < Y_dim; qq++){
          prob_temp(q) = prob_temp(q) -0.5*log(sigma(q, qq)) - 0.5*(Y(h, qq) - mu(q, qq))*(Y(h, qq) - mu(q, qq))/sigma(q, qq) ;
        }
        for (qq = 0; qq < V_dim; qq++){
          prob_temp(q) = prob_temp(q) + weight(h)*V(h, qq)*log(omega(q, qq) + 0.00000000000000000001);
        }
        sum_temp = sum_temp + prob_temp(q); 
      }
      
      for(q = 0; q < K; q++)
      {
        prob_temp(q) = prob_temp(q) - sum_temp/K; 
        if (prob_temp(q) > 709){prob_temp(q) = 709;}
        if (prob_temp(q) < -300){prob_temp(q) = -300; }
        prob_temp(q) = exp(prob_temp(q));
        sum_temp2 = sum_temp2 + prob_temp(q);
      }
      
      if(sum_temp2 == 0){
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = 1.0/K;
        }
        Rcout <<i<< "% zero\n";
      }
      else{
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = prob_temp(q)/sum_temp2;
        }
      }
      
      Z(h) = sample(sample_index, 1, 1==0, prob_temp)(0); 
    }
    
    
    // Monitor the process
    if (i*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        mu_store(i, count) = mu(q, qq);
        sigma_store(i, count) = sigma(q, qq);
        count++;
      }
    }
    
    count = 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < V_dim; qq++)
      {
        omega_store(i, count) = omega(q, qq);
        count++;
      }
    }
    
    for(q = 0; q < N; q++)
    {
      Z_store(i, q) = Z(q);
    }
    
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu_store, Rcpp::Named("sigma") = sigma_store, Rcpp::Named("omega") = omega_store, Rcpp::Named("Z") = Z_store);
}



// [[Rcpp::export]]
arma::mat rdirichlet_cpp(NumericVector alpha, int size) {

  arma::mat sample = arma::zeros(1, size);
  
  double sum_temp = 0;
  for (int j = 0; j < size; ++j){
    double cur = rgamma(1, alpha(j), 1.0)(0);
    sample(0, j) = cur; 
    sum_temp += cur; 
  }
  
  for (int j = 0; j < size; ++j){
    sample(0, j) = sample(0, j)/sum_temp; 
  }
  return(sample);
}

// [[Rcpp::export]]
double lbeta_function(arma::mat alpha) {
  double re = 0.0;
  double sum_temp = 0.0;
  
  int N =alpha.n_cols;
  for (int i = 0; i < N; i++){
    re += lgamma_function(exp(alpha(0, i)));
    sum_temp += exp(alpha(0, i));
  }
  
  re -= lgamma_function(sum_temp);

  return(re);
  }
  
  
// [[Rcpp::export]]
double lgamma_function(double alpha) {
  NumericVector alpha_vec(1);
  
  alpha_vec(0) = alpha;
  double qq = Rcpp::sum(Rcpp::lgamma(alpha_vec));
  
  return(qq);
}


// [[Rcpp::export]]
double dnorm_function(arma::mat alpha, NumericVector mu, arma::mat sigma) {
  arma::vec qq = dmvnorm(alpha, mu, sigma);
  
  return(qq(0, 0));
}





// [[Rcpp::export]]
arma::mat rnorm_function(int n, NumericVector mu, arma::mat sigma) {
  arma::mat qq = rmvnorm(n, mu, sigma);
  
  return(qq);
}






