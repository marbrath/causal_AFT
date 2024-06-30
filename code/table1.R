library("extRemes")
library("survival")
library("statmod")
library("gamma")

datgen_u0 <- function(num, rho0, HR0, u0_dist) {
  epsilon = mu + sigma*(-revd(num, loc = 0, scale = 1)) #note -1 bc. parametrization of evd in R

  if (u0_dist=="gamma") {
    u0 = rgamma(num, shape=1/rho0, scale=rho0)
  }

  if (u0_dist=="invgauss") {
    u0 = rinvgauss(num, mean=1, disp=rho0)
  }

  vv_a = -log(u0)*sigma - log(HR0)*sigma + epsilon
  vv_0 = -log(u0)*sigma + epsilon
  a = rbinom(num, 1, 0.5)
  vv = -log(u0)*sigma - log(HR0)*sigma*a + epsilon
  tt = exp(vv)

  df = data.frame(time = tt, log_time_a = vv_a, log_time_0 = vv_0, status = 1, log_time = vv, a = a)
  return(df)
}

lambda_ = 1/60 # Weibull scale
k = 3 # Weibull shape

mu = -log(lambda_)/k
sigma = 1/k

rho0_values = c(0.5, 1, 2)
HR0_values = c(1/3., 3)

num_rho0 = length(rho0_values)
num_HR0 = length(HR0_values)

u0_dist = "gamma"

n_sim = 10000
n_obs = 500

cox_res = array(numeric(), c(num_rho0, num_HR0, n_sim))
osd_res = array(numeric(), c(num_rho0, num_HR0, n_sim))

res_table_cox = array(numeric(), c(num_rho0, num_HR0))
res_table_osd = array(numeric(), c(num_rho0, num_HR0))
rmsd_table_cox = array(numeric(), c(num_rho0, num_HR0))
rmsd_table_osd = array(numeric(), c(num_rho0, num_HR0))

ci_table_cox = array(numeric(),c(num_rho0, num_HR0, 2))
ci_table_osd = array(numeric(),c(num_rho0, num_HR0, 2))


set.seed(1)

for (rho0_idx in 1:num_rho0)
{
  rho0 = rho0_values[rho0_idx];

  for (HR0_idx in 1:num_HR0)
  {
    HR0 = HR0_values[HR0_idx];

    for (i in 1:n_sim) {
      print(i)
      df = datgen_u0(n_obs, rho0, HR0, u0_dist)
      
      s_cox = coxph(Surv(time, status)~a, data=df)
      
      osd_res[rho0_idx,HR0_idx,i] = mean(df$log_time[df$a==0]) - mean(df$log_time[df$a==1])
      cox_res[rho0_idx,HR0_idx,i] = s_cox$coefficients[1]
    }
  }
}

for (rho0_idx in 1:num_rho0) {
  for (HR0_idx in 1:num_HR0) {
    theta = log(HR0_values[HR0_idx])*sigma
    print(theta)
    
    res_table_cox[rho0_idx,HR0_idx] = round(exp(mean(cox_res[rho0_idx,HR0_idx,])),3)
    res_table_osd[rho0_idx,HR0_idx] = round(exp(mean(osd_res[rho0_idx,HR0_idx,])),3)
     
    rmsd_table_cox[rho0_idx,HR0_idx] = sqrt(mean((cox_res[rho0_idx,HR0_idx,] - log(HR0_values[HR0_idx]))**2))/sqrt(n_sim)
    rmsd_table_osd[rho0_idx,HR0_idx] = sqrt(mean((osd_res[rho0_idx,HR0_idx,] - theta)**2))/sqrt(n_sim)
    
    l_cox = exp(mean(cox_res[rho0_idx,HR0_idx,]) - 2*rmsd_table_cox[rho0_idx,HR0_idx])
    u_cox = exp(mean(cox_res[rho0_idx,HR0_idx,]) + 2*rmsd_table_cox[rho0_idx,HR0_idx])
    l_osd = exp(mean(osd_res[rho0_idx,HR0_idx,]) - 2*rmsd_table_osd[rho0_idx,HR0_idx])
    u_osd = exp(mean(osd_res[rho0_idx,HR0_idx,]) + 2*rmsd_table_osd[rho0_idx,HR0_idx])
    
    ci_table_cox[rho0_idx,HR0_idx,1] = round(l_cox,3)
    ci_table_cox[rho0_idx,HR0_idx,2] = round(u_cox,3)
    ci_table_osd[rho0_idx,HR0_idx,1] = round(l_osd,3)
    ci_table_osd[rho0_idx,HR0_idx,2] = round(u_osd,3)
  }
}

