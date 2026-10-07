// ---- Syphilis in women of child-bearing age, gestational syphilis and CS ----
// Sensitivity analysis (delay distribution): model_primary.stan with the calendar-period
// structure of the diagnosis-to-delivery delay distribution passed as data: delay_period[t]
// in 1..P assigns quarter t to one of P periods, each with its own simplex[3] p_delay[p].
// Run with P=1 (a single distribution) and P=2 (2016-2021 vs. 2022-2025). The primary model
// corresponds to P=3 (2016-2019, 2020-2021, 2022-2025).
// Expects q_base (K,T) as baseline fertility per woman,
// rr[c] scales fertility by age-category.
//
// The CS risk r(t) is supplied as data (p_r_fixed): p_r[1] for 2016-2021 and p_r[2] for
// 2022-2025. The primary analysis uses 11/76 and 28/375; the CS-risk sensitivity analysis
// re-runs the model at 10 Latin hypercube points drawn from Beta(12, 66) and Beta(29, 348).
// The lower bound of the fertility-rate ratio s(a) (rr) is also supplied as data (rr_floor;
// 1 in the primary analysis).
//
// Delay from diagnosis to delivery: three-point distribution, estimated separately for
// 2016-2019, 2020-2021 and 2022-2025.
// Scenario analysis: cum_cs_w[w] gives cumulative CS when a proportion scen_mult[w] of women
// with syphilis has the fertility rate of all women in the same age group.

data {
  int<lower=1> T;                  // 2016Q1 ... (e.g., to 2025Q3 -> T=39)
  int<lower=1> C;                  // # age categories (e.g., 6)
  int<lower=1> K;                  // # age indices (fine grid)
  int<lower=1> Y;                  // # year (2016-2025)

  array[C] int<lower=1, upper=K> lower_bound;  // category -> k-range lower
  array[C] int<lower=1, upper=K> upper_bound;  // category -> k-range upper
  array[K] int<lower=1, upper=C> cat_of_k;     // each k belongs to which category c

  array[C, T] int Case;     // reported cases by category/time
  array[T]    int CS;      // CS counts by time

  // Externalized fertility baselines (already divided by 1000):
  array[K, T]   real<lower=0> q_base;     // for 1..T

  // Externalized yearly births (already divided by 100000):
  array[Y] real<lower=0> brt_per100000;

  // scenario assessment
  int<lower=1> W;                       // # of scenarios
  array[W] real<lower=0, upper=1> scen_mult;  // proportion in each scenario

  // CS risk among viable deliveries: p_r[1] applies to 2016-2021 (t=1..24),
  // p_r[2] applies to 2022-2025 (t=25..T).
  vector<lower=0,upper=1>[2] p_r_fixed;

  // Lower bound of rr (fertility-rate ratio s(a)): 1 in the primary analysis; 0.8 and 0.5
  // (and 0 with a lognormal prior, model_rrprior.stan) in the sensitivity analyses.
  real<lower=0, upper=1> rr_floor;

  int<lower=1> P;                                  // number of delay-distribution periods
  array[T] int<lower=1, upper=P> delay_period;     // period of each quarter
}

parameters {
  array[T] real<lower=0> i_t;     // time trend of incidence
  real<lower=0> s_i;              // RW step sd on log(i_t)
  simplex[K] i_a;                 // age mixing 
  array[P] simplex[3] p_delay;   // delay probability, one simplex per period
  // rr >= rr_floor (= 1 in the primary analysis; see Methods for the structural reason)
  vector<lower=rr_floor>[C] rr;   // age-cat fertility multipliers (replaces rr1..rr6)
  real<lower=0, upper=1> rep;              // reporting rate 
}

transformed parameters {
  // pmf of diagnosis-to-report delay
  // g[delay k, time t] 
array[K, T] real g;
for (t in 1:T) {
  for (k in 1:K) {
    if (k <= 3) g[k, t] = p_delay[delay_period[t], k];
    else        g[k, t] = 0.0;
  }
}

  // CS risk supplied as data
  vector[2] p_r = p_r_fixed;

  // time-varying CS risk among viable deliveries
  array[T] real r;
  for (t in 1:24) r[t] = p_r[1];   // 2016-2021
  for (t in 25:T) r[t] = p_r[2];   // 2022-2025
  
  // expected notified cases by age/time from diagnosed infection
  array[K, T] real e_case;
  for (t in 1:T) {
    for (k in 1:K) {
      e_case[k, t] = i_t[t] * i_a[k] * rep;
    }
  }
  
  // aggregate to age-categories
  array[C, T] real e_case_c;
  for (t in 1:T) {
    for (c in 1:C) {
      e_case_c[c, t] = sum(e_case[lower_bound[c]:upper_bound[c], t]);
    }
  }
  
  // CS 
  array[K, T] real e_cs_a;//expected notified CS by age/time from diagnosed infection convolved with reporting delay g
  for (t in 1:T) {
    for (k in 1:K) {
      real sum_temp = 0;
      int max_s = (k > t) ? t : k;
      for (s in 1:max_s) {
        int kk = k - (s - 1);
        int tt = t - (s - 1);
        int c = cat_of_k[kk];
        if (kk >= 1 && tt >= 1) sum_temp += i_t[tt] * i_a[kk] * q_base[kk, tt]* rr[c] * r[tt] * g[s,tt];
      }
      e_cs_a[k, t] = sum_temp;
    }
  }

  array[T] real e_cs; //expected notified CS by time
  for (t in 1:T) e_cs[t] = sum(e_cs_a[1:K, t]);
  
  
}

model {
  // likelihoods
  for (t in 1:T)
    for (c in 1:C)
      Case[c, t] ~ poisson(e_case_c[c, t]);
  
  for (t in 1:T)
    CS[t] ~ poisson(e_cs[t]);
  
  // RW prior on i_t
  i_t[1] ~ lognormal(log(10), 1);
  for (t in 2:T)
    log(i_t[t]) ~ normal(log(i_t[t-1]), s_i);
  
  i_a ~ dirichlet(rep_vector(1.0, K));
 
  s_i ~ normal(0, 1); // half-normal due to lower=0 constraint

  // reporting fraction: no prior beyond its [0, 1] bounds

  for (p in 1:P)
    p_delay[p] ~ dirichlet(rep_vector(1.0, 3));

  // rr: no prior beyond its lower bound rr_floor
}

generated quantities {
  // ratio of the 2022-2025 to the 2016-2021 CS risk
  real rrisk = p_r[2] / p_r[1];

  // diagnosed cases related
  array[C, T] real i_c; // by age category/time
  for (t in 1:T)
    for (c in 1:C)
      i_c[c, t] = i_t[t] * sum(i_a[lower_bound[c]:upper_bound[c]]);
  
  // CS related 
  array[Y] real cs_pop; //aggregate to year per 100000 births
  for (y in 1:Y){
    if (y <= 9) cs_pop[y] = brt_per100000[y]*sum(e_cs[(4*y-3):(4*y)]);
    else cs_pop[y] = brt_per100000[y]*sum(e_cs[(4*y-3):(4*y-1)]);
  }

  
  // gestational syphilis related
  array[K, T] real ip_a;// incident gestational syphilis by age/time
  for (t in 1:T) {
    for (k in 1:K) {
      int c = cat_of_k[k];
      ip_a[k, t] = i_t[t] * i_a[k] * q_base[k, t] * rr[c];
    }
  }
  
  array[C, T] real ip_c;// aggregate to age-categories
  for (t in 1:T) {
    for (c in 1:C) {
      ip_c[c, t] = sum(ip_a[lower_bound[c]:upper_bound[c], t]);
    }
  }
  
  array[C, Y] real gst_y;// aggregate to year by age-categories
  for (c in 1:C){
    for (y in 1:Y){
    if (y <= 9) gst_y[c, y] = sum(ip_c[c, (4*y-3):(4*y)]);
    else gst_y[c, y] = sum(ip_c[c, (4*y-3):(4*y-1)]);
  }
  }
  
  array[T] real ip_t;// aggregate to time 
  for (t in 1:T) {
      ip_t[t] = sum(ip_c[1:C, t]);
  }
  
  array[Y] real gst_pop;// aggregate to year per 100000 births
  for (y in 1:Y){
    if (y <= 9) gst_pop[y] = brt_per100000[y]*sum(ip_t[(4*y-3):(4*y)]);
    else gst_pop[y] = brt_per100000[y]*sum(ip_t[(4*y-3):(4*y-1)]);
  }
  
  // posterior predicted notified cases by category/time
  array[C, T] int case_prd;
  for (t in 1:T)
    for (c in 1:C)
      case_prd[c, t] = poisson_rng(e_case_c[c, t]);
   
  // posterior predicted CS by time   
  array[T] int cs_prd; 
  for (t in 1:T) cs_prd[t] = poisson_rng(e_cs[t]);
  
  // simulated cumulative CS
  array[W, K, T] real e_cs_a_w;//expected notified CS by age/time from diagnosed infection convolved with reporting delay g
  for (w in 1:W){
    for (t in 1:T) {
    for (k in 1:K) {
      real sum_temp = 0;
      int max_s = (k > t) ? t : k;
      for (s in 1:max_s) {
        int kk = k - (s - 1);
        int tt = t - (s - 1);
        int c = cat_of_k[kk];
        if (kk >= 1 && tt >= 1) sum_temp += i_t[tt] * i_a[kk] * q_base[kk, tt]* ((1-scen_mult[w])*rr[c]+scen_mult[w]*1) * r[tt] * g[s,tt];
      }
      e_cs_a_w[w, k, t] = sum_temp;
    }
  }}
  
  array[W,T] real e_cs_w; //expected notified CS by time
  for (w in 1:W){
    for (t in 1:T) e_cs_w[w,t] = sum(e_cs_a_w[w, 1:K, t]);
  }
  
  array[W,T] int cs_prd_w; // posterior predicted CS by time 
  for (w in 1:W){
  for (t in 1:T) cs_prd_w[w,t] = poisson_rng(e_cs_w[w,t]);
  }
  
  array[W] real cum_cs_w; // posterior predicted CS from 2016-2025
  for (w in 1:W){
    cum_cs_w[w] = sum(cs_prd_w[w, 1:T]);
  }
  
  
  // log-lik for WAIC/LOO
  array[C, T] real log_lik_case;
  for (c in 1:C)
    for (t in 1:T)
      log_lik_case[c, t] = poisson_lpmf(Case[c, t] | e_case_c[c, t]);
  
  vector[T] log_lik_cs;
  for (t in 1:T)
    log_lik_cs[t] = poisson_lpmf(CS[t] | e_cs[t]);
  
}
