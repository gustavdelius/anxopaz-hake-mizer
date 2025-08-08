
#include <TMB.hpp>

template<class Type>
vector<Type> calculate_F_mort(Type logit_l50, Type log_ratio_left, Type log_l50_right_offset, Type log_ratio_right,
                              Type log_catchability, vector<Type> blength, Type min_len, Type max_len) 
{
  Type l50 = min_len + (max_len - min_len) * invlogit(logit_l50);  
  Type l25 = l50 * (1 - exp(log_ratio_left)); 
  
  Type l50_right = l50 + exp(log_l50_right_offset); 
  Type l25_right = l50_right * (1 + exp(log_ratio_right));
  
  Type catchability = exp(log_catchability);
  
  vector<Type> F_mort(blength.size());
  
  Type sr = l50 - l25;
  Type s1 = l50 * log(Type(3.0)) / sr;
  Type s2 = s1 / l50;
  
  Type sr_right = l50_right - l25_right;
  Type s1_right = l50_right * log(Type(3.0)) / sr_right;
  Type s2_right = s1_right / l50_right;
  
  for (int i = 0; i < blength.size(); i++) {
    Type ilength = blength(i);
    Type sel = (Type(1.0) / (Type(1.0) + exp(s1 - s2 * ilength))) * 
      (Type(1.0) / (Type(1.0) + exp(s1_right - s2_right * ilength)));
    F_mort(i) = catchability * sel;
  }
  
  return F_mort;
}

template<class Type>
vector<Type> calculate_natmort(vector<Type> weight, Type log_M, Type d)
{
  Type M = exp(log_M);
  vector<Type> natmort = M * pow(weight, d);
  return natmort;
}



template<class Type>
vector<Type> calculate_mort(vector<Type> total_F_mort, vector<Type> natmort, vector<Type> weight)
{
  vector<Type> mort = natmort + total_F_mort;
  return mort;
}


template<class Type>
vector<Type> calculate_erepro(vector<Type> weight, Type log_h, Type log_n, Type log_ks, Type log_p, Type log_k, 
                              Type log_f0, Type log_alpha)
{
  Type h = exp(log_h);
  Type n = Type(1.0) / (Type(1.0) + exp(-log_n));
  Type ks = exp(log_ks);
  Type p = exp(log_p);
  Type k = exp(log_k);
  Type alpha = Type(1.0) / (Type(1.0) + exp(-log_alpha));
  Type f0 = Type(1.0) / (Type(1.0) + exp(-log_f0));
  vector<Type> encounter = h * pow(weight, n);
  vector<Type> metab = ks * pow(weight, p) + k * weight;
  vector<Type> erepro = alpha * f0 * encounter - metab;
  return erepro;
}


template<class Type>
vector<Type> calculate_repro_prop(vector<Type> weight, Type log_w_mat, Type log_n)
{
  Type w_mat = exp(log_w_mat);
  Type n = Type(1.0) / (Type(1.0) + exp(-log_n));
  vector<Type> reprop = pow(weight / w_mat, Type(1.0) - n);
  
  Type max_val = reprop[0];
  for(int i = 1; i < reprop.size(); i++) {
    if(reprop[i] > max_val) max_val = reprop[i];
  }
  
  vector<Type> repro_prop = reprop / max_val;
  
  return repro_prop;
}


template<class Type>
vector<Type> calculate_growth(vector<Type> erepro, vector<Type> repro_prop, Type log_w_mat, Type log_U, vector<Type> weight)
{
  Type w_mat = exp(log_w_mat);
  Type U = exp(log_U);
  Type c1 = Type(1.0);
  vector<Type> psi = repro_prop / (c1 + pow(weight / w_mat, -U));
  vector<Type> growth = erepro * (c1 - psi);
  return growth;
}


template<class Type>
vector<Type> calculate_N(vector<Type> mort, vector<Type> growth, Type biomass, vector<Type> bin_widths,
                         vector<Type> weight)
{
  int size = bin_widths.size();
  vector<Type> N(size);
  N(0) = Type(1.0);
  for (int i = 1; i < size; ++i) {
    Type denominator = growth(i) + mort(i) * bin_widths(i);
    N(i) = N(i - 1) * growth(i - 1) / denominator;
  }
  Type total_biomass = Type(0.0);
  for (int i = 0; i < size; ++i) {
    total_biomass += N(i) * bin_widths(i) * weight(i);
  }
  N = N * biomass / total_biomass;
  
  return N;
}

template<class Type>
vector<Type> calculate_catch_per_bin(vector<Type> N, vector<Type> F_mort, vector<Type> bin_widths)
{
  
  vector<Type> densities = N * F_mort;
  int num_bins = bin_widths.size();
  vector<Type> catch_per_bin(num_bins);
  for (int i = 0; i < num_bins; ++i) {
    catch_per_bin[i] = bin_widths[i] * densities[i];
  }
  return catch_per_bin;
}

template<class Type>
Type calculate_yield(vector<Type> catch_per_bin, vector<Type> weight)
{
  Type model_yield = Type(0.0);
  for (int i = 0; i < catch_per_bin.size(); ++i) {
    model_yield += catch_per_bin[i] * weight[i];
  }
  return model_yield;
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  // **Data Section**
  DATA_MATRIX(counts);
  DATA_VECTOR(bin_widths);
  DATA_VECTOR(weight);
  DATA_VECTOR(blength);
  DATA_VECTOR(yield);
  DATA_SCALAR(biomass);
  DATA_SCALAR(min_len);
  DATA_SCALAR(max_len);
  DATA_SCALAR(yield_lambda);
  DATA_SCALAR(biomass_lambda);
  DATA_INTEGER(n_g);
  
  // **Parameter Section**
  PARAMETER_VECTOR(logit_l50);
  PARAMETER_VECTOR(log_ratio_left);
  PARAMETER_VECTOR(log_l50_right_offset);
  PARAMETER_VECTOR(log_ratio_right);
  PARAMETER_VECTOR(log_catchability); 
  
  PARAMETER(log_U);
  PARAMETER(log_M);
  PARAMETER(d);
  PARAMETER(log_h);
  PARAMETER(log_n);
  PARAMETER(log_ks);
  PARAMETER(log_p);
  PARAMETER(log_k);
  PARAMETER(log_f0);
  PARAMETER(log_alpha);
  PARAMETER(log_w_mat);
  
  
  int n_bins = bin_widths.size();
  
  vector<Type> total_F_mort(n_bins);
  total_F_mort.setZero(); 
  matrix<Type> F_mort_mat(n_bins, n_g);
  for (int g = 0; g < n_g; ++g) {
    vector<Type> F_mort_g = calculate_F_mort(
      logit_l50[g], log_ratio_left[g], log_l50_right_offset[g], log_ratio_right[g], log_catchability[g], blength, min_len, max_len);
    for (int i = 0; i < n_bins; ++i) {
      F_mort_mat(i, g) = F_mort_g(i);
      total_F_mort(i) += F_mort_g(i);
    }
  }
  
  vector<Type> mort = calculate_mort(total_F_mort, log_M, d, weight);
  
  vector<Type> erepro = calculate_erepro(weight, log_h, log_n, log_ks, log_p, log_k, log_f0, log_alpha);
  
  vector<Type> repro_prop = calculate_repro_prop(weight, log_w_mat, log_n);
  
  vector<Type> growth = calculate_growth(erepro, repro_prop, log_w_mat, log_U, weight);
  
  vector<Type> N = calculate_N(mort, growth, biomass, bin_widths, weight);
  
  matrix<Type> catch_per_bin_mat(n_bins, n_g);
  
  for (int g = 0; g < n_g; ++g) {
    
    Eigen::Matrix<Type, Eigen::Dynamic, 1> F_mort_column = F_mort_mat.col(g);
    
    vector<Type> F_mort_g(n_bins);
    for (int i = 0; i < n_bins; i++) {
      F_mort_g(i) = F_mort_column(i);
    }
    
    vector<Type> catch_per_bin_g = calculate_catch_per_bin(N, F_mort_g, bin_widths);
    catch_per_bin_mat.col(g) = catch_per_bin_g;
    
  }
  
  vector<Type> model_yield_g(n_g);
  
  for (int g = 0; g < n_g; ++g) {
    
    Eigen::Matrix<Type, Eigen::Dynamic, 1> catch_per_bin_column = catch_per_bin_mat.col(g);
    
    vector<Type> catch_per_bin_g(n_bins);
    for (int i = 0; i < n_bins; i++) {
      catch_per_bin_g(i) = catch_per_bin_column(i);
    }
    
    Type yield_g = calculate_yield(catch_per_bin_g, weight);
    model_yield_g(g) = yield_g;
    
  }
  
  Type nll = Type(0.0);
  
  for (int g = 0; g < n_g; ++g) {
    vector<Type> catch_per_bin_g = catch_per_bin_mat.col(g);
    vector<Type> counts_g = counts.col(g);
    
    vector<Type> probs_g = catch_per_bin_g + Type(1e-10);
    probs_g = probs_g / probs_g.sum();
    
    Type nll_g = -dmultinom(counts_g, probs_g, true);
    
    nll_g += yield_lambda * pow(log(model_yield_g(g) / yield(g)), Type(2));
    
    nll += nll_g;
  }
  
  
  // Type biomass_model = Type(0.0);
  // for (int i = 0; i < N.size(); ++i) {
  //   biomass_model += N(i) * bin_widths(i) * weight(i);
  // }
  // nll += biomass_lambda * pow(log(biomass_model / biomass), Type(2));
  
  // TMBAD_ASSERT(nll >= 0);
  TMBAD_ASSERT(CppAD::isfinite(nll));
  if (!CppAD::isfinite(nll)) error("nll is not finite");
  
  return nll;
}
