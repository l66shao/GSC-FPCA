

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Function to sample from a multivariate normal distribution
arma::vec rmvnorm(const arma::vec& mu, const arma::mat& Sigma) {
  int n = mu.n_elem;
  arma::vec z = arma::randn(n);
  return mu + arma::chol(Sigma, "lower") * z;
}

// [[Rcpp::export]] 

List run_mcmc(int nloop, int burnin, int thin, List Y_sparse, List time_sparse,  List X_sparse,List X_sparse_c,
              List phi_t, arma::vec Mu,arma::mat Beta,arma::vec Beta_c, arma::mat Alpha, arma::mat Theta, 
              arma::mat D, arma::vec Sigma_omega, arma::mat Sigma_alpha, 
              double c_omega, int est_mu, int est_alpha, int est_d, 
              int est_sigma_omega, int est_theta,int est_beta,int est_beta_c,int PSLOO,arma::vec b_control,int loo_treated_mean=2) {
  
  double eta1_eps=0.01;
  double eta2_eps=0.01;
  
  int N = Y_sparse.size();
  int N_c = arma::sum(b_control);
  int P = 1;
  int sum_K = Mu.n_elem;  
  int sum_Q = D.n_rows;//number of PCs
  
  arma::mat loglike(N,(nloop-burnin),arma::fill::zeros);
  
  List ALPHA_array, MU_array, D_array, SIGMA_OMEGA_array, THETA_array,BETA_array,BETA_c_array,ALPHA_TM_array,NU_array;
  arma::vec Nu = arma::zeros(sum_Q);
  arma::mat Alpha_tm = Alpha;
  int array_ind = 0;
  
  for (int iter = 0; iter < nloop; iter++) {
    
    // ----------------Update Mu---------------
    if (est_mu == 1) {
      arma::mat Sigma_mu = (1.0 / c_omega) * arma::eye(sum_K, sum_K);
      arma::vec mu_mu = arma::zeros(sum_K);
      
      int cols = Beta.n_cols;
      for (int i = 0; i < N; i++) {
        if(b_control[i]==1){
        arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
        arma::mat X_i= as<arma::mat>(X_sparse[i]);
        arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
        arma::mat Sigma_omega_i = Sigma_omega(0) * arma::eye(Y_i.n_elem, Y_i.n_elem);
        arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
        
        arma::vec tmp = arma::zeros(Y_i.n_elem);  // Initialize tmp to zero vector of same length as Y_i
        
        for (int ii = 0; ii < cols; ++ii) {
          arma::mat Dt = diagmat(X_i.col(ii));                   // D is (n × n), diagonal matrix
          arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);     // Result: (n × 1)
          tmp += contrib;                                       // Accumulate across covariates
        }
        arma::vec Y_i_vec_resid = Y_i - Phi_D_t_i * Theta * Alpha.col(i);
        arma::mat inv_term = arma::inv(Phi_D_t_i * Theta * Sigma_alpha * trans(Theta) * trans(Phi_D_t_i) + Sigma_omega_i);
        Sigma_mu += trans(Phi_D_t_i) * inv_term * Phi_D_t_i;
        mu_mu += trans(Phi_D_t_i) * inv_term * (Y_i-tmp-X_c_i*Beta_c);
        
     }
      }
      
      Sigma_mu = 0.5 * (Sigma_mu + Sigma_mu.t());
      Sigma_mu = arma::inv(Sigma_mu);
      mu_mu = Sigma_mu * mu_mu;
      Mu = rmvnorm(mu_mu, Sigma_mu);
    }
    
    // ----------------Update Beta----------------
    
    if (est_beta == 1){
      int cols = Beta.n_cols;
      for (int jj = 0; jj < cols; ++jj) {
      arma::mat Sigma_beta = (1.0 / c_omega) * arma::eye(Beta.n_rows, Beta.n_rows);
      arma::vec mu_beta = arma::zeros(Beta.n_rows);
     
     
          for (int i = 0; i < N; i++) {
    if(b_control[i]==1){
            arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
            arma::mat X_i= as<arma::mat>(X_sparse[i]);
            arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
            arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
            arma::mat Sigma_omega_i = Sigma_omega(0) * arma::eye(Y_i.n_elem, Y_i.n_elem);
            
            arma::vec Y_i_vec_resid = Y_i - Phi_D_t_i * Theta * Alpha.col(i);
            //10*1-10*7* 7*2 *2*1
            arma::mat inv_term = arma::inv(Phi_D_t_i * Theta * Sigma_alpha * trans(Theta) * trans(Phi_D_t_i) + Sigma_omega_i);
            //need to check this M here
            arma::vec tmp = arma::zeros(Y_i.n_elem);  // initialize tmp with zeros
            for (int ii = 0; ii < cols; ++ii) {
              if (ii == jj) continue;  // skip diagonal
              arma::mat Dt = diagmat(X_i.col(ii));            // D: (n × n)
              arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);  // result: (n × 1)
              tmp += contrib;
            }
            
            Sigma_beta +=trans(diagmat(X_i.col(jj))*(Phi_D_t_i)) * inv_term * (diagmat(X_i.col(jj))*(Phi_D_t_i));
            //2*10 * 10*10 * 10*2
           
            mu_beta +=trans(diagmat(X_i.col(jj))*Phi_D_t_i)  * inv_term * (Y_i-(Phi_D_t_i)*Mu-tmp-X_c_i*Beta_c);//here need to add
           }
          }
          
      Sigma_beta = 0.5 * (Sigma_beta + Sigma_beta.t());
      Sigma_beta = arma::inv(Sigma_beta);
      mu_beta = Sigma_beta * mu_beta;
      Beta.col(jj) = rmvnorm(mu_beta, Sigma_beta);
      }
    }
    
    // ----------------Update Beta_c----------------
    if (est_beta_c == 1) {
      
      arma::mat Sigma_beta_c = (1.0 / c_omega) * arma::eye(Beta_c.size(), Beta_c.size());
      arma::vec mu_beta_c = arma::zeros(Beta_c.size());
      
      for (int i = 0; i < N; i++) {
        if(b_control[i]!=1) continue;
        arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
        arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
        arma::mat X_i= as<arma::mat>(X_sparse[i]);
        
        arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
        arma::mat Sigma_omega_i = Sigma_omega(0) * arma::eye(Y_i.n_elem, Y_i.n_elem);
        
        arma::vec Y_i_vec_resid = Y_i - Phi_D_t_i * Theta * Alpha.col(i);
        //10*1-10*7* 7*2 *2*1
        arma::mat inv_term = arma::inv(Phi_D_t_i * Theta * Sigma_alpha * trans(Theta) * trans(Phi_D_t_i) + Sigma_omega_i);
        
        Sigma_beta_c += trans(X_c_i) * inv_term * X_c_i;
        //3*10 * 10*10 * 10*2, T_i=10
        int cols = Beta.n_cols;
        arma::vec tmp = arma::zeros(Y_i.n_elem);  // initialize tmp with zeros
        for (int ii = 0; ii < cols; ++ii) {
          arma::mat Dt = diagmat(X_i.col(ii));            // D: (n × n)
          arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);  // result: (n × 1)
          tmp += contrib;
        }
        
        mu_beta_c += trans(X_c_i) * inv_term * (Y_i-(Phi_D_t_i)*Mu-tmp);
        
      }
      
      Sigma_beta_c = 0.5 * (Sigma_beta_c + Sigma_beta_c.t());
      Sigma_beta_c = arma::inv(Sigma_beta_c);
      mu_beta_c = Sigma_beta_c * mu_beta_c;
      Beta_c = rmvnorm(mu_beta_c, Sigma_beta_c);
      
    }
    
    
    
    // ----------------Update Alpha----------------
    if (est_alpha == 1) {
      int cols = Beta.n_cols;
      for (int i = 0; i < N; i++) {
        if(b_control[i]!=1) continue;  // treated scores are drawn after the Theta update, below
        
        arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
        arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
        arma::mat Sigma_omega_i = Sigma_omega(0) * arma::eye(Y_i.n_elem, Y_i.n_elem);
        arma::mat X_i= as<arma::mat>(X_sparse[i]);
        
        arma::vec tmp = arma::zeros(Y_i.n_elem);  // Initialize tmp to zero vector of same length as Y_i
        arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
        
        for (int ii = 0; ii < cols; ++ii) {
          arma::mat Dt = diagmat(X_i.col(ii));                   // Dt is (n × n), diagonal matrix
          arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);     // Result: (n × 1)
          tmp += contrib;                                       // Accumulate across covariates
        }
        arma::vec Y_i_vec_resid = - (Y_i - Phi_D_t_i * Mu  - tmp - X_c_i*Beta_c);
        arma::mat Sigma_alpha_i = trans(Theta) * trans(Phi_D_t_i) * arma::inv(Sigma_omega_i) * Phi_D_t_i * Theta + arma::inv(Sigma_alpha);
        Sigma_alpha_i = 0.5 * (Sigma_alpha_i + Sigma_alpha_i.t());
        Sigma_alpha_i = arma::inv(Sigma_alpha_i);
        arma::vec mu_alpha_i = trans(Theta) * trans(Phi_D_t_i) * arma::inv(Sigma_omega_i) * -Y_i_vec_resid;
        mu_alpha_i = Sigma_alpha_i * mu_alpha_i;
        Alpha.col(i) = rmvnorm(mu_alpha_i, Sigma_alpha_i);
        
        }
    }
    
    // ----------------Update D----------------
    if (est_d == 1) {
      int cols = Beta.n_cols;
     arma::mat alpha_tmp=Alpha.cols(arma::find(b_control == 1));
      for (int j = 0; j < sum_Q; j++) {
        double eta1_D_j_p=eta1_eps;
        double eta2_D_j_p=eta2_eps;
         eta1_D_j_p = 0.5 * N_c + eta1_D_j_p;
         eta2_D_j_p = eta2_D_j_p + arma::accu(arma::square(alpha_tmp.row(j))) / 2.0;
        D(j, j) = 1.0 / R::rgamma(eta1_D_j_p, 1.0 / eta2_D_j_p);
      }
    }
    
    // ----------------Update Sigma_alpha----------------
    Sigma_alpha = D;
    
    // ----------------Update Sigma_omega----------------
    if (est_sigma_omega == 1) {
      int indk = 0, indq = 0;
      int cols = Beta.n_cols;
      for (int p = 0; p < P; p++) {
        double eta1_omega_p = eta1_eps;
        double eta2_omega_p =eta2_eps;
        for (int i = 0; i < N; i++) {
          if(b_control[i]==1){
          arma::vec Y_pi = as<arma::vec>(Y_sparse[i]);
          eta1_omega_p +=  0.5*Y_pi.n_elem;
          arma::mat Phi_D_p = as<arma::mat>(phi_t[i]);
          arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
          arma::mat X_i= as<arma::mat>(X_sparse[i]);
          arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
          arma::vec tmp = arma::zeros(Y_pi.n_elem);  // Initialize tmp to zero vector of same length as Y_i
          
          for (int ii = 0; ii < cols; ++ii) {
            arma::mat Dt = diagmat(X_i.col(ii));                   // Dt is (n × n), diagonal matrix
            arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);     // Result: (n × 1)
            tmp += contrib;                                       // Accumulate across covariates
          }
          
          arma::vec Y_pi_resid  = Y_pi -tmp- Phi_D_t_i * Theta * Alpha.col(i)- Phi_D_p.t()* Mu -X_c_i*Beta_c;;
          
          
          
          //arma::vec Y_pi_resid=Phi_D_p.t();
          eta2_omega_p += 0.5 * arma::accu(arma::square(Y_pi_resid));
          }}
        Sigma_omega(p) = 1.0 / R::rgamma(eta1_omega_p, 1.0 / eta2_omega_p);
        indk += sum_K / P;
        indq += sum_Q / P;
      }
    }
    
    // ----------------Update Theta----------------
    if (est_theta == 1) {
      int cols = Beta.n_cols;
      int indk = 0, indq = 0;
      for (int p = 0; p < P; p++) {
        arma::vec mu_theta_p = arma::zeros(sum_Q * sum_K);//2*7
        arma::mat Sigma_theta_p = (1.0 / c_omega) * arma::eye(sum_Q * sum_K, sum_Q * sum_K);//14*14
        
        for (int i = 0; i < N; i++) {
          
          if(b_control[i]==1){
          arma::vec Y_pi = as<arma::vec>(Y_sparse[i]);
          arma::mat X_i= as<arma::mat>(X_sparse[i]);
          arma::vec Alpha_pi = Alpha.col(i);
          arma::mat Alpha_pim = arma::diagmat(Alpha.col(i));
         
         arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
         
         
          for (int t = 0; t < Y_pi.n_elem; t++) {//t stands for time
            arma::mat Phi_D_pit =arma::zeros(sum_Q*sum_K, sum_Q);
            arma::vec tmp = arma::zeros(Y_pi.n_cols);  // Initialize tmp to zero vector of same length as Y_i
         
            for( int q=0; q<sum_Q;q++){
              Phi_D_pit.submat(q*sum_K, q, (q+1)*sum_K-1, q) = as<arma::mat>(phi_t[i]).col(t);
            }
            for (int ii = 0; ii < cols; ++ii) {// Dt is (n × n), diagonal matrix
              arma::vec contrib = X_i(t,ii) * trans(as<arma::mat>(phi_t[i]).col(t)) * Beta.col(ii);     // Result: (n × 1)
              tmp += contrib;                                       // Accumulate across covariates
            }
            
            arma::vec Y_pit_resid = Y_pi.row(t) -tmp- trans(as<arma::mat>(phi_t[i]).col(t)) * Mu- X_c_i.row(t)*Beta_c;
            Sigma_theta_p += Phi_D_pit * Alpha_pi* trans(Alpha_pi) * trans(Phi_D_pit) / Sigma_omega(p);
            mu_theta_p += Phi_D_pit * Alpha_pi * Y_pit_resid / Sigma_omega(p);
          }
          }}
        
        Sigma_theta_p = arma::inv(0.5 * (Sigma_theta_p + Sigma_theta_p.t()));
        Sigma_theta_p=0.5 * (Sigma_theta_p + Sigma_theta_p.t());
        mu_theta_p = Sigma_theta_p * mu_theta_p;
        Theta.submat(indk, indq, sum_K / P - 1, sum_Q / P - 1) = reshape(rmvnorm(mu_theta_p, Sigma_theta_p), sum_K / P, sum_Q / P);
        
        indk += sum_K / P;
        indq += sum_Q / P;
      }
    }
    
    // ----------------Treated scores and treated-group mean----------------
    // The components above come from controls only. Treated units enter here and nowhere else.
    // nu is the treated group's mean score, with a flat prior and the individual scores integrated out:
    //   nu | components, pre-exposure data ~ N(G^{-1} b, G^{-1}), G = sum_i D^{-1}(I - R_i), b = sum_i D^{-1} a_i,
    // where a_i and Omega_i are the posterior mean and covariance of alpha_i under prior N(0, D), and
    // R_i = Omega_i D^{-1}. Under prior N(nu, D) the posterior mean is a_i + R_i nu with the same Omega_i.
    // Alpha holds draws under prior mean 0 and Alpha_tm draws under prior mean nu, with common noise.
    if (est_alpha == 1) {
      int cols = Beta.n_cols;
      arma::mat Dinv = arma::inv(Sigma_alpha);
      arma::mat G = arma::zeros(sum_Q, sum_Q);
      arma::vec b = arma::zeros(sum_Q);
      std::vector<arma::mat> R_list;
      for (int i = 0; i < N; i++) {
        if(b_control[i]==1) continue;
        arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
        arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
        arma::mat X_i= as<arma::mat>(X_sparse[i]);
        arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
        arma::vec tmp = arma::zeros(Y_i.n_elem);
        for (int ii = 0; ii < cols; ++ii) {
          tmp += diagmat(X_i.col(ii)) * Phi_D_t_i * Beta.col(ii);
        }
        arma::vec u_i = Y_i - Phi_D_t_i * Mu - tmp - X_c_i * Beta_c;
        arma::mat F_i = Phi_D_t_i * Theta;
        arma::mat Prec_i = trans(F_i) * F_i / Sigma_omega(0) + Dinv;
        Prec_i = 0.5 * (Prec_i + Prec_i.t());
        arma::mat Omega_i = arma::inv(Prec_i);
        Omega_i = 0.5 * (Omega_i + Omega_i.t());
        arma::vec a_i = Omega_i * trans(F_i) * u_i / Sigma_omega(0);
        arma::mat R_i = Omega_i * Dinv;
        G += Dinv * (arma::eye(sum_Q, sum_Q) - R_i);
        b += Dinv * a_i;
        R_list.push_back(R_i);
        Alpha.col(i) = rmvnorm(a_i, Omega_i);
      }
      if (R_list.size() > 0) {
        G = 0.5 * (G + G.t());
        arma::mat G_inv = arma::inv(G);
        G_inv = 0.5 * (G_inv + G_inv.t());
        Nu = rmvnorm(G_inv * b, G_inv);
        Alpha_tm = Alpha;
        int ind_tr = 0;
        for (int i = 0; i < N; i++) {
          if(b_control[i]==1) continue;
          Alpha_tm.col(i) = Alpha.col(i) + R_list[ind_tr] * Nu;
          ind_tr++;
        }
      }
    }

    // Leave-one-out log predictive densities.
    // Controls: log N(u_i; 0, V_i), with the scores integrated out, for PSIS.
    // Treated units (loo_treated_mean == 2): their data enter only through nu, so the leave-one-out predictive density
    // given the components is exact: u_i | components, other treated data ~ N(F_i nu_{-i}, V_i + F_i G_{-i}^{-1} F_i'), with
    // G_{-i} = G - F_i' V_i^{-1} F_i and nu_{-i} = G_{-i}^{-1} (b - F_i' V_i^{-1} u_i). The components do not depend on the
    // treated data, so averaging this density over the draws needs no importance weights.
    // loo_treated_mean == 1 gives the density given the current nu, and 0 the density under prior mean zero.
    if (PSLOO == 1&&iter % thin == 0 && iter >= burnin) {
      arma::mat G_loo = arma::zeros(sum_Q, sum_Q);
      arma::vec b_loo = arma::zeros(sum_Q);
      std::vector<int> tr_i;
      std::vector<arma::mat> tr_G, tr_F, tr_V;
      std::vector<arma::vec> tr_b, tr_u;
      arma::mat Sigma_mu = (1.0 / c_omega) * arma::eye(sum_K, sum_K);
      arma::vec mu_mu = arma::zeros(sum_K);
      
      int cols = Beta.n_cols;
      for (int i = 0; i < N; i++) {
        //if(b_control[i]==1){
        arma::vec Y_i = as<arma::vec>(Y_sparse[i]);
        arma::mat X_i= as<arma::mat>(X_sparse[i]);
        arma::mat X_c_i= as<arma::mat>(X_sparse_c[i]);
        arma::mat Phi_D_t_i = trans(as<arma::mat>(phi_t[i]));
        arma::mat Sigma_omega_i = Sigma_omega(0) * arma::eye(Y_i.n_elem, Y_i.n_elem);
        
        arma::vec tmp = arma::zeros(Y_i.n_elem);  // Initialize tmp to zero vector of same length as Y_i
        
        for (int ii = 0; ii < cols; ++ii) {
          arma::mat Dt = diagmat(X_i.col(ii));                   // Dt is (n × n), diagonal matrix
          arma::vec contrib = Dt * Phi_D_t_i * Beta.col(ii);     // Result: (n × 1)
          tmp += contrib;                                       // Accumulate across covariates
        }
        //- Phi_D_t_i * Theta * Alpha.col(i)
        arma::vec Y_i_vec_resid = Y_i - Phi_D_t_i * Mu  - tmp-X_c_i*Beta_c;
        if(b_control[i]!=1 && loo_treated_mean==1) Y_i_vec_resid -= Phi_D_t_i * Theta * Nu;
        arma::mat inv_term = (Phi_D_t_i * Theta * Sigma_alpha * trans(Theta) * trans(Phi_D_t_i) + Sigma_omega_i);
        if(b_control[i]!=1 && loo_treated_mean==2){
          arma::mat F_i = Phi_D_t_i * Theta;
          arma::mat Vinv_F = arma::solve(inv_term, F_i);
          arma::mat G_i = trans(F_i) * Vinv_F;
          arma::vec b_i = trans(Vinv_F) * Y_i_vec_resid;
          G_loo += G_i; b_loo += b_i;
          tr_i.push_back(i); tr_G.push_back(G_i); tr_F.push_back(F_i); tr_V.push_back(inv_term); tr_b.push_back(b_i); tr_u.push_back(Y_i_vec_resid);
          continue;
        }
        double likei=arma::as_scalar(exp(-0.5*trans(Y_i_vec_resid)*arma::inv(inv_term)*Y_i_vec_resid));
        likei=likei/sqrt(arma::det(inv_term));
        likei=likei/(std::pow(2*M_PI,Y_i.n_elem/2.0));
        loglike(i,(iter-burnin)) = log(likei) ;
        
        // }
      }
      for (size_t j = 0; j < tr_i.size(); j++) {
        arma::mat G_mi = G_loo - tr_G[j];
        G_mi = 0.5 * (G_mi + G_mi.t());
        arma::mat G_mi_inv = arma::inv(G_mi);
        arma::vec nu_mi = G_mi_inv * (b_loo - tr_b[j]);
        arma::vec r = tr_u[j] - tr_F[j] * nu_mi;
        arma::mat S = tr_V[j] + tr_F[j] * G_mi_inv * trans(tr_F[j]);
        S = 0.5 * (S + S.t());
        double ld, sg;
        arma::log_det(ld, sg, S);
        loglike(tr_i[j],(iter-burnin)) = -0.5 * (r.n_elem * log(2 * M_PI) + ld + arma::as_scalar(trans(r) * arma::solve(S, r)));
      }
      
    }
    
    
    // Store results after burn-in
    if (iter % thin == 0 && iter >= burnin) {
      ALPHA_array.push_back(Alpha);
      MU_array.push_back(Mu);
      D_array.push_back(D);
      SIGMA_OMEGA_array.push_back(Sigma_omega);
      THETA_array.push_back(Theta);
      BETA_array.push_back(Beta);
      BETA_c_array.push_back(Beta_c);
      ALPHA_TM_array.push_back(Alpha_tm);
      NU_array.push_back(Nu);
    }
  }
  
  return List::create(ALPHA_array, MU_array, D_array, D_array, SIGMA_OMEGA_array, THETA_array,BETA_array,BETA_c_array,loglike,ALPHA_TM_array,NU_array);
}
