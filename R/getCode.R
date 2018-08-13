
getCode <- function(modelName){

  model.Po <- paste0('data{
  int J; // number of animals
  int ystararaw[J]; // after treatment McMaster count
  int ystarbraw[J]; // before treatment McMaster count
  int fpre[J];
  int fpost[J];
  real w[J]; // weights
  real wmo;  // weighted mean of the outliers
  real postmean; // mean of the post egg counts
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
  real<lower=0> mub[J];
  real<lower=1> alpha;
}
transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  for (i in 1:J){
   lambdab[i] = mub[i]/fpre[i];
   lambdaa[i] = delta*mub[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);   // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  alpha ~ normal(wmo/postmean,10);
  mub ~ gamma(kappa, kappa/mu);  // likelihoods
  ystarbraw ~ poisson(lambdab);
  for (n in 1:J){
   target += w[n]*poisson_lpmf(ystararaw[n] | lambdaa[n]) +
          (1-w[n])*poisson_lpmf(ystararaw[n] | alpha*lambdaa[n] );
}
}
')

  model.ZIPo <- paste0('data{
  int J; // number of animals
  int ystararaw[J]; // after treatment McMaster count
  int ystarbraw[J]; // before treatment McMaster count
  int fpre[J];
  int fpost[J];
  real w[J]; // weights
  real wmo;  // weighted mean of the outliers
  real postmean; // mean of the post egg counts
  }
  parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
  real<lower=0> mub[J];
  real<lower=0,upper=1> phi;
  real<lower=1> alpha;
  }
  transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  for (i in 1:J){
  lambdab[i] = mub[i]/fpre[i];
  lambdaa[i] = delta*mub[i]/fpost[i];
  }
  }
  model{
  mu ~ gamma(1,0.001);   // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  phi ~ beta(1,1);
  alpha ~ normal(wmo/postmean,10);
  mub ~ gamma(kappa, kappa/mu);        // likelihoods
  for (n in 1:J){
  if (ystarbraw[n] == 0)
  target += log_sum_exp(bernoulli_lpmf(1 | phi),
  bernoulli_lpmf(0 | phi)+poisson_lpmf(ystarbraw[n] | lambdab[n]));
  else
  target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarbraw[n] | lambdab[n]);
  }
  for (n in 1:J){
  if (ystararaw[n] == 0)
  target += log_sum_exp(bernoulli_lpmf(1 | phi),
  bernoulli_lpmf(0 | phi)+poisson_lpmf(ystararaw[n] | lambdaa[n]));
  else
  target += bernoulli_lpmf(0 | phi) +  w[n]*poisson_lpmf(ystararaw[n] | lambdaa[n])+
  (1-w[n])*poisson_lpmf(ystararaw[n] | alpha*lambdaa[n] );
  }
  }
  ')

  model.UPo <- paste0('data{
  int Ja; // number of animals
  int Jb;
  int ystararaw[Ja]; // after treatment McMaster count
  int ystarbraw[Jb]; // before treatment McMaster count
  int fpre[Ja];
  int fpost[Jb];
  real w[Ja]; // weights
  real wmo;  // weighted mean of the outliers
  real postmean; // mean of the post egg counts
  }
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0, upper = 1> delta;
  real<lower=0> mub[Jb]; # true epg before treatment
  real<lower=0> mua[Ja]; # true epg after treatment
  real<lower=1> alpha;
}
transformed parameters{
  real lambdaa[Ja];
  real lambdab[Jb];
  for (i in 1:Jb){
   lambdab[i] = mub[i]/fpre[i];
  }
  for (i in 1:Ja){
   lambdaa[i] = delta*mua[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);   // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  alpha ~ normal(wmo/postmean,10);
  mub ~ gamma(kappa, kappa/mu)+
  mua ~ gamma(kappa, kappa/mu);  // likelihoods
  ystarbraw ~ poisson(lambdab);
  for (n in 1:Ja){
   target += w[n]*poisson_lpmf(ystararaw[n] | lambdaa[n]) +
          (1-w[n])*poisson_lpmf(ystararaw[n] | alpha*lambdaa[n] );
  }
}
  ')

  model.ZIUPo <- paste0('data{
  int Ja; // number of animals
  int Jb;
  int ystararaw[Ja]; // after treatment McMaster count
  int ystarbraw[Jb]; // before treatment McMaster count
  int fpre[Ja];
  int fpost[Jb];
  real w[Ja]; // weights
  real wmo;  // weighted mean of the outliers
  real postmean; // mean of the post egg counts
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0, upper=1> delta;
  real<lower=0> mub[Jb];
  real<lower=0> mua[Ja];
  real<lower=0,upper=1> phi;
  real<lower=1> alpha;
}
transformed parameters{
  real lambdaa[Ja];
  real lambdab[Jb];
  for (i in 1:Jb){
    lambdab[i] = mub[i]/fpre[i];
  }
  for (i in 1:Ja){
    lambdaa[i] = delta*mua[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);   // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  phi ~ beta(1,1);
  alpha ~ normal(wmo/postmean,10);
  mub ~ gamma(kappa, kappa/mu);  // likelihoods
  mua ~ gamma(kappa, kappa/mu);
  for (n in 1:Jb){
    if (ystarbraw[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | phi),
                bernoulli_lpmf(0 | phi)+poisson_lpmf(ystarbraw[n] | lambdab[n]));
    else
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarbraw[n] | lambdab[n]);
  }
  for (n in 1:Ja){
    if (ystararaw[n] == 0){
      target += log_sum_exp(bernoulli_lpmf(1 | phi),
    bernoulli_lpmf(0 | phi)+poisson_lpmf(ystararaw[n] | lambdaa[n]));
    }else{
    target += bernoulli_lpmf(0 | phi) + w[n]*poisson_lpmf(ystararaw[n] | lambdaa[n]) +
              (1-w[n])*poisson_lpmf(ystararaw[n]| alpha*lambdaa[n]);
    }
  }
}
  ')

  switch (modelName,
   "Po" = model.Po,
   "ZIPo" = model.ZIPo,
   "UPo" = model.UPo,
   "ZIUPo" = model.ZIUPo
  )
}
