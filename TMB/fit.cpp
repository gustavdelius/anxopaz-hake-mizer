
#include <TMB.hpp>

template<class Type> vector<Type> wlf(vector<Type> w, Type log_a, Type log_b){
  Type a = exp(log_a);
  Type b = exp(log_b);
  vector<Type> L = pow(w/a, Type(1.0)/b);
  return L;
}

template<class Type> vector<Type> lwf(vector<Type> l, Type log_a, Type log_b){
  Type a = exp(log_a);
  Type b = exp(log_b);
  vector<Type> W = a * pow(l, b);
  return W;
}

template<class Type> vector<Type> diff_w(vector<Type> w) {
  int n = w.size();
  vector<Type> dw(n);
  for(int i = 0; i < n - 1; i++) {
    dw(i) = w(i+1) - w(i);
  }
  dw(n-1) = Type(0);
  return dw;
}

template<class Type> vector<Type> Nres(vector<Type> wp, Type kappa, Type lambda){
  vector<Type> Nr(wp.size());
  for(int i=0; i<wp.size(); i++){
    Nr(i) = kappa * pow(wp(i), -lambda);
  }
  return Nr;
}


template<class Type> matrix<Type> pred_kernel(vector<Type> w, vector<Type> wp, Type log_beta, Type log_sigma){
  Type beta = exp(log_beta);
  Type sigma = exp(log_sigma);
  int nw = w.size();
  int nwp = wp.size();
  matrix<Type> predk(nw, nwp);
  for(int i=0; i<nw; i++){
    for(int j=0; j<nwp; j++){
      predk(i,j) = exp(-pow(log(w(i)/wp(j)) - log(beta), 2) / (Type(2.0) * sigma*sigma));
    }
  }
  return predk;
}


template<class Type>
vector<Type> Eavailable(vector<Type> w, vector<Type> wp, vector<Type> dw, vector<Type> dwp,
                Type inter_HR, Type inter_HH, vector<Type> N, vector<Type> Nr, matrix<Type> predk){
  int nw = w.size();
  vector<Type> Eavail(nw);
  for(int i=0; i<nw; i++){
    Type enc_re = 0.0;
    Type enc_os = 0.0;
    for(int j=0; j<wp.size(); j++){
      enc_re += predk(i,j) * inter_HR * Nr(j) * wp(j) * dwp(j);
    }
    for(int j=0; j<w.size(); j++){
      enc_os += predk(j,i) * inter_HH * N(j) * w(j) * dw(j);
    }
    Eavail(i) = enc_re + enc_os;
  }
  return Eavail;
}

template<class Type>
vector<Type> Eencounter(vector<Type> w, vector<Type> Eavail, Type log_gamma, Type log_q){
  Type q = exp(log_q);
  Type gamma = exp(log_gamma);
  vector<Type> Eenc(w.size());
  for(int i=0; i<w.size(); i++){
    Eenc(i) = Eavail(i) * gamma * pow(w(i), q);
  }
  return Eenc;
}

template<class Type>
vector<Type> feeding_level(vector<Type> w, vector<Type> Eenc, Type log_h, Type log_n){
  Type h = exp(log_h);
  Type n = Type(1.0) / (Type(1.0) + exp(-log_n));
  vector<Type> feed_level(w.size());
  for(int i=0; i<w.size(); i++){
    feed_level(i) = Eenc(i) / (Eenc(i) + h * pow(w(i), n));
  }
  return feed_level;
}
  

template<class Type>
vector<Type> Ereproandgrowth(vector<Type> w, Type log_h, Type log_n, Type log_ks, Type log_p, Type log_k,
                             Type log_alpha, vector<Type> feed_level){
  Type h = exp(log_h);
  Type n = Type(1.0) / (Type(1.0) + exp(-log_n));
  Type ks = exp(log_ks);
  Type p = exp(log_p);
  Type k = exp(log_k);
  Type alpha = Type(1.0) / (Type(1.0) + exp(-log_alpha));
  
  vector<Type> erg(w.size());
  for(int i=0; i<w.size(); i++){
    Type emet = ks * pow(w(i), p);
    Type emov = k * w(i);
    Type maxc = h * pow(w(i), n);
    erg(i) = alpha * feed_level(i) * maxc - emet - emov;
    if(erg(i) < 0) erg(i) = 0;
  }
  return erg;
}

template<class Type>
vector<Type> calculate_psi(vector<Type> w, Type log_w_max, Type log_w_mat, Type log_U, Type log_n){
  
  Type U = exp(log_U);
  Type w_mat = exp(log_w_mat);
  Type w_max = exp(log_w_max);
  Type n = Type(1.0) / (Type(1.0) + exp(-log_n));
  
  vector<Type> psi(w.size());
  for(int i=0; i<w.size(); i++){
    Type repp = pow(w(i)/w_max, 1.0 - n);
    Type mat = 1.0 / (1.0 + pow(w(i)/w_mat, -U));
    psi(i) = repp * mat;
  }
  return psi;
}

template<class Type>
vector<Type> PredM(vector<Type> w, vector<Type> dw, Type inter_HH, Type log_gamma, Type log_q,
                   vector<Type> N, vector<Type> feed_level, matrix<Type> predk){
  Type q = exp(log_q);
  Type gamma = exp(log_gamma);
  vector<Type> pmort(w.size());
  for(int i=0; i<w.size(); i++){
    pmort(i) = 0.0;
    for(int j=0; j<w.size(); j++){
      pmort(i) += inter_HH * predk(j,i) * (1.0 - feed_level(j)) * gamma * pow(w(j), q) * N(j) * dw(j);
    }
  }
  return pmort;
}


template<class Type>
vector<Type> calculate_F_mort(Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
    Type log_catchability, Type log_effort, vector<Type> L, Type min_len, Type max_len) 
{
  Type l50 = min_len + (max_len - min_len) * invlogit(logit_l50);  
  Type l25 = l50 * (1 - exp(log_ratio_left)); 
  
  Type l50_right = l50 + exp(log_l50_right_offset); 
  Type l25_right = l50_right * (1 + exp(log_ratio_right));
  
  Type catchability = exp(log_catchability);
  Type effr = exp(log_effort);
  
  vector<Type> F_mort(L.size());
  
  Type sr = l50 - l25;
  Type s1 = l50 * log(Type(3.0)) / sr;
  Type s2 = s1 / l50;
  
  Type sr_right = l50_right - l25_right;
  Type s1_right = l50_right * log(Type(3.0)) / sr_right;
  Type s2_right = s1_right / l50_right;
  
  for (int i = 0; i < L.size(); i++) {
    Type ilength = L(i);
    Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * ilength))) * 
      (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * ilength)));
    F_mort(i) = catchability * sel * effr;
  }
  
  return F_mort;
}


template<class Type>
vector<Type> nat_mort(vector<Type> w, Type log_M, Type d){
  Type M = exp(log_M);
  vector<Type> nmort = M * pow(w, d);
  return nmort;
}


template<class Type>
vector<Type> calculate_mort(vector<Type> total_F_mort, vector<Type> pmort, vector<Type> nmort){
  vector<Type> mort = pmort + nmort + total_F_mort;
  return mort;
}


template<class Type>
vector<Type> Egrowth(vector<Type> erg, vector<Type> psi){
  
  vector<Type> growth = erg * (Type(1.0) - psi);
  return growth;
}


template<class Type>
vector<Type> calculate_N(vector<Type> w, vector<Type> dw, vector<Type> mort, vector<Type> growth, 
                         Type biomass)
{
  int size = dw.size();
  vector<Type> N(size);
  N(0) = Type(1.0);
  for (int i = 1; i < size; ++i) {
    Type denominator = growth(i) + mort(i) * dw(i);
    N(i) = N(i - 1) * growth(i - 1) / denominator;
  }
  Type total_biomass = Type(0.0);
  for (int i = 0; i < size; ++i) {
    total_biomass += N(i) * dw(i) * w(i);
  }
  N = N * biomass / total_biomass;
  
  return N;
}

template<class Type>
vector<Type> calculate_catch_per_bin(vector<Type> N, vector<Type> F_mort, vector<Type> dw)
{
  
  vector<Type> densities = N * F_mort;
  int num_bins = dw.size();
  vector<Type> catch_per_bin(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    catch_per_bin[i] = dw[i] * densities[i];
  }
  return catch_per_bin;
}

template<class Type>
Type calculate_yield(vector<Type> catch_per_bin, vector<Type> w)
{
  Type model_yield = Type(0.0);
  for (int i = 0; i < catch_per_bin.size(); ++i) {
    model_yield += catch_per_bin[i] * w[i];
  }
  return model_yield;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // **Data Section**
  DATA_MATRIX(counts);
  DATA_VECTOR(w);
  DATA_VECTOR(wp);
  DATA_VECTOR(yield);
  DATA_SCALAR(biomass);
  
  DATA_SCALAR(yield_lambda);
  DATA_SCALAR(biomass_lambda);
  DATA_INTEGER(n_g);
  
  // **Parameter Section**
  PARAMETER_VECTOR(logit_l50);
  PARAMETER_VECTOR(log_ratio_left);
  PARAMETER_VECTOR(log_l50_right_offset);
  PARAMETER_VECTOR(log_ratio_right);
  PARAMETER_VECTOR(log_catchability); 
  PARAMETER_VECTOR(log_effort); 
  
  PARAMETER(kappa);
  PARAMETER(lambda);
  PARAMETER(log_a);
  PARAMETER(log_b);
  PARAMETER(log_beta);
  PARAMETER(log_sigma);
  PARAMETER(inter_HR);
  PARAMETER(inter_HH);
  PARAMETER(log_gamma);
  PARAMETER(log_q);
  PARAMETER(log_h);
  PARAMETER(log_n);
  PARAMETER(log_ks);
  PARAMETER(log_p);
  PARAMETER(log_k);
  PARAMETER(log_alpha);
  PARAMETER(log_U);
  PARAMETER(log_w_mat);
  PARAMETER(log_w_min);
  PARAMETER(log_w_max);
  PARAMETER(log_M);
  PARAMETER(d);
  
  int n_bins = w.size();
  
  vector<Type> dw = diff_w(w);
  vector<Type> dwp = diff_w(wp);
  
  vector<Type> Nr = Nres(wp, kappa, lambda);
  matrix<Type> predk = pred_kernel(w, wp, log_beta, log_sigma);
  
  matrix<Type> F_mort_mat(n_bins, n_g);
  vector<Type> total_F_mort(n_bins);
  
  vector<Type> N(n_bins);

  Type w_min = exp(log_w_min);
  Type w_max = exp(log_w_max);

  Type min_len = pow(w_min/exp(log_a), Type(1.0)/exp(log_b));
  Type max_len = pow(w_max/exp(log_a), Type(1.0)/exp(log_b));
  vector<Type> L = wlf( w, log_a, log_b);
  
  for(int i=0; i<n_bins; i++) {
    N(i) = biomass / (n_bins * w(i)); 
  }
  
  int max_iter = 50;
  Type tolerance = 1e-6;
  vector<Type> N_old(n_bins);
  
  for(int iter = 0; iter < max_iter; iter++) {
    N_old = N;
    
    vector<Type> Eavail = Eavailable(w, wp, dw, dwp, inter_HR, inter_HH, N, Nr, predk);
    vector<Type> Eenc = Eencounter(w, Eavail, log_gamma, log_q);
    vector<Type> feed_level = feeding_level(w, Eenc, log_h, log_n);
    vector<Type> erg = Ereproandgrowth(w, log_h, log_n, log_ks, log_p, log_k, log_alpha, feed_level);
    vector<Type> psi = calculate_psi(w, log_w_max, log_w_mat, log_U, log_n);
    vector<Type> pmort = PredM(w, dw, inter_HH, log_gamma, log_q, N, feed_level, predk);
    
    total_F_mort.setZero();
    for (int g = 0; g < n_g; ++g) {
      vector<Type> F_mort_g = calculate_F_mort(logit_l50[g], log_ratio_left[g], log_l50_right_offset[g], 
           log_ratio_right[g], log_catchability[g], log_effort[g], L, min_len, max_len);
      for (int i = 0; i < n_bins; ++i) {
        F_mort_mat(i, g) = F_mort_g(i);
        total_F_mort(i) += F_mort_g(i);
      }
    }
    
    vector<Type> nmort = nat_mort(w, log_M, d);
    vector<Type> mort = calculate_mort(total_F_mort, pmort, nmort);
    vector<Type> growth = Egrowth(erg, psi);
    
    N = calculate_N(w, dw, mort, growth, biomass);
    
    Type diff = sqrt((N - N_old).square().sum() / n_bins);
    if(diff < tolerance) {
      break;
    }
  }

  vector<Type> Eavail = Eavailable(w, wp, dw, dwp, inter_HR, inter_HH, N, Nr, predk);
    vector<Type> Eenc = Eencounter(w, Eavail, log_gamma, log_q);
    vector<Type> feed_level = feeding_level(w, Eenc, log_h, log_n);
    vector<Type> erg = Ereproandgrowth(w, log_h, log_n, log_ks, log_p, log_k, log_alpha, feed_level);
    vector<Type> psi = calculate_psi(w, log_w_max, log_w_mat, log_U, log_n);
    vector<Type> pmort = PredM(w, dw, inter_HH, log_gamma, log_q, N, feed_level, predk);
    
    total_F_mort.setZero();
    for (int g = 0; g < n_g; ++g) {
      vector<Type> F_mort_g = calculate_F_mort(logit_l50[g], log_ratio_left[g], log_l50_right_offset[g], 
           log_ratio_right[g], log_catchability[g], log_effort[g], L, min_len, max_len);
      for (int i = 0; i < n_bins; ++i) {
        F_mort_mat(i, g) = F_mort_g(i);
        total_F_mort(i) += F_mort_g(i);
      }
    }
    
    vector<Type> nmort = nat_mort(w, log_M, d);
    vector<Type> mort = calculate_mort(total_F_mort, pmort, nmort);
    vector<Type> growth = Egrowth(erg, psi);
    
    N = calculate_N(w, dw, mort, growth, biomass);
  
  matrix<Type> catch_per_bin_mat(n_bins, n_g);
  vector<Type> model_yield_g(n_g);
  
  for (int g = 0; g < n_g; ++g) {
    vector<Type> F_mort_g(n_bins);
    for (int i = 0; i < n_bins; i++) {
      F_mort_g(i) = F_mort_mat(i, g); 
    }
    
    vector<Type> catch_per_bin_g = calculate_catch_per_bin(N, F_mort_g, dw);
    catch_per_bin_mat.col(g) = catch_per_bin_g;
    model_yield_g(g) = calculate_yield(catch_per_bin_g, w);
  }
  
  Type nll = Type(0.0);
  
  for (int g = 0; g < n_g; ++g) {
    vector<Type> catch_per_bin_g = catch_per_bin_mat.col(g);
    vector<Type> counts_g = counts.col(g);
    
    vector<Type> probs_g = catch_per_bin_g + Type(1e-10);
    probs_g = probs_g / probs_g.sum();
    nll -= dmultinom(counts_g, probs_g, true);
    
    nll += yield_lambda * pow(log(model_yield_g(g) / yield(g)), Type(2));
  }
  
  Type biomass_model = Type(0.0);
  for (int i = 0; i < N.size(); ++i) {
    biomass_model += N(i) * dw(i) * w(i);
  }
  nll += biomass_lambda * pow(log(biomass_model / biomass), Type(2));
 
  TMBAD_ASSERT(CppAD::isfinite(nll));
  if (!CppAD::isfinite(nll)) error("nll is not finite");
  
  vector<Type> final_B(n_bins);
  for (int i = 0; i < N.size(); ++i) {
    final_B(i) = N(i) * dw(i) * w(i);
  }
  
  REPORT(final_B);
  return nll;
}
