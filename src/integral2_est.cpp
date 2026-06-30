#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;


// [[Rcpp::export]]
SEXP integral2_est(
    arma::vec init_parameters,
    arma::vec time,
    arma::vec event,
    arma::vec IV,
    arma::vec IV_c,
    arma::mat Covariates,
    arma::mat D_status,
    arma::vec stime,
    int max_iter,
    double tol,
    double eta,
    double contraction
)
{
  int n = time.size();
  int k = stime.size();
  int p = init_parameters.size() - 1;
  // printf("%d %d %d", n, k, p);
  
  double betaD, Asyvar, hatsigma, hatphi, ave_hatsigma;
  int step = 0;
  bool Convergence;
  arma::vec hatsigma_v(n), beta(p), res(n), Covbeta(n), fn(p+1), SY(k), dSY(k),A(k), B(k), dLam_dbetaD(k), dLam(k), pk(p+1), new_parameters(p+1), diff(p+1), fn_abs(p+1);
  arma::mat int_D(n, k), dNt(n, k), Yt(n, k), K(n, k), int_expbetaD(n, k), dLam_dbeta(p, k), int_dexpbetaD(n, k), dPhi(n, p+1), Hessian(p+1, p+1);
  // double H_com_val, bktrkg;
  
  for (int iter = 0; iter < max_iter; iter++)
  {
    
    res.zeros(); Covbeta.zeros(), dPhi.zeros(), Hessian.zeros();
    fn.zeros(); SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(), int_dexpbetaD.zeros();
    A.zeros(), B.zeros();
    betaD = init_parameters[0];
    beta = init_parameters.subvec(1, p); // index is only adapted to the one dimensional covariates
    
    // arma::arma_print(beta);
    // if(iter == 0) 
    // {
    //     arma::arma_print(betaD);
    //     arma::arma_print(beta);
    // }
    for (int i = 0; i < n; i++)
    {
      for (int ii = 0; ii < p; ii++)
      {
        Covbeta[i] += Covariates(i, ii) * beta[ii];
      }
    }
    
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
      {
        if(j == 0)
          int_D(i, j) = IV[i] * stime[j];
        else 
          int_D(i, j) = int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]);
        
        
        dNt(i, j) = (time[i] == stime[j]) ? 1:0;
        dNt(i, j) *= event[i];
        Yt(i, j) = (time[i] >= stime[j]) ? 1:0;
        SY[j] = SY[j] + exp(betaD * int_D(i, j)) * Yt(i, j);
        dSY[j] = dSY[j] + int_D(i, j) * exp(betaD * int_D(i, j)) * Yt(i, j);
      }
      
    }
    // arma::arma_print(SY);
    // printf("dd\n");
    
    dLam.zeros();
    dLam_dbetaD.zeros();
    dLam_dbeta.zeros();
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
      {
        if (j == 0)
        {
          dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*stime[j] - Covbeta[i]*stime[j])/ SY[j];
          dLam_dbetaD[j] = dLam_dbetaD[j]+Yt(i, j)*int_D(i, j)*exp(int_D(i, j)*betaD)*(dNt(i,j) - (D_status(i,j)*betaD + Covbeta[i])*stime[j])/SY[j] - 
            Yt(i, j)*exp(int_D(i, j)*betaD)*D_status(i,j)*stime[j]/SY[j] - 
            dSY[j]*Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*stime[j] - Covbeta[i]*stime[j])/ (SY[j] * SY[j]);
          for (int kk = 0; kk < p; kk++)
          {
            dLam_dbeta(kk, j) = dLam_dbeta(kk, j) - Yt(i, j)*exp(int_D(i, j)*betaD)*Covariates(i, kk)*stime[j]/SY[j];
          }
        } else
        {
          dLam[j] = dLam[j] + Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*(stime[j] - stime[j-1]) -
            Covbeta[i]*(stime[j] - stime[j-1]))/ SY[j];
          dLam_dbetaD[j] = dLam_dbetaD[j]+Yt(i, j)*int_D(i, j)*exp(int_D(i, j)*betaD)*(dNt(i,j) - (D_status(i,j)*betaD + Covbeta[i])*(stime[j] - stime[j-1]))/SY[j] - 
            Yt(i, j)*exp(int_D(i, j)*betaD)*D_status(i,j)*(stime[j] - stime[j-1])/SY[j] - 
            dSY[j]*Yt(i, j) * exp(int_D(i, j) * betaD) * (dNt(i, j) - D_status(i, j)*betaD*(stime[j] - stime[j-1])- Covbeta[i]*(stime[j] - stime[j-1]))/ (SY[j] * SY[j]);
          for (int kk = 0; kk < p; kk++)
          {
            dLam_dbeta(kk, j) = dLam_dbeta(kk, j) - Yt(i, j)*exp(int_D(i, j)*betaD)*Covariates(i, kk)*(stime[j] - stime[j-1])/SY[j];
          }
        }
      }
    }
    // printf("2\n");
    
    
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k; j++)
      {
        // Integral computation for exp(D\theta) and D * exp(D\theta)
        if (j == 0)
        {
          if (betaD != 0)
          {
            int_expbetaD(i, j) = (IV[i] > 0.5 ) ? (exp(betaD * int_D(i, j))-1)/betaD:stime[j];
            int_dexpbetaD(i, j) = (IV[i] > 0.5) ? (int_D(i, j) * exp(betaD * int_D(i, j))/betaD - int_expbetaD(i, j)/betaD):0;
          }
          else
          { 
            int_expbetaD(i, j) = stime[j];
            int_dexpbetaD(i, j) = int_D(i, j) * int_D(i, j) / 2;
          }
          res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j)- Yt(i, j) * dLam[j]) - int_expbetaD(i, j) * (IV[i] * betaD + Covbeta[i]);
          
          for (int kk = 0; kk < p+1; kk++)
          {
            if (kk == 0)
            {
              
              dPhi(i, kk) = dPhi(i, kk) + (dNt(i, j) - Yt(i, j) * dLam[j]) * int_D(i, j) * exp(betaD * int_D(i, j)) -
                exp(betaD*int_D(i,j))*Yt(i,j)*dLam_dbetaD[j]-
                int_dexpbetaD(i, j) * (IV[j] * betaD + Covbeta[i]) -
                int_expbetaD(i, j) * IV[j];
            }
            else
            {
              dPhi(i, kk) = dPhi(i, kk) - Covariates(i, kk-1) * int_expbetaD(i, j)-
                Yt(i, j) * exp(int_D(i, j)*betaD)*dLam_dbeta(kk-1, j);
            }
            
          }
        }
        else
        {
          if (betaD != 0)
          {
            int_expbetaD(i, j) = (D_status(i, j-1) > 0.5) ? (exp(betaD * int_D(i, j)) - exp(betaD * int_D(i, j-1)))/betaD:(stime[j] - stime[j-1]) * exp(betaD * int_D(i, j));
            int_dexpbetaD(i, j) = (D_status(i, j-1) > 0.5) ? ((int_D(i, j) * exp(betaD * int_D(i, j)) - int_D(i, j-1) * exp(betaD * int_D(i, j-1)))/betaD - int_expbetaD(i, j)/betaD) : (int_D(i, j)*exp(betaD * int_D(i, j))*(stime[j] - stime[j-1]));
          }
          else 
          {
            int_expbetaD(i, j) = stime[j] - stime[j-1];
            int_dexpbetaD(i, j) = (stime[j] - stime[j-1]) * int_D(i, j-1) + D_status(i, j-1) * (stime[j] - stime[j-1]) * (stime[j] - stime[j-1])/2;
            
          }
          res[i] = res[i] + exp(int_D(i, j) * betaD) * (dNt(i, j) - Yt(i, j)*dLam[j]) - Yt(i, j-1) * int_expbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]);
          
          for (int kk = 0; kk < p+1; kk++)
          {
            if (kk == 0)
            {
              
              dPhi(i, kk) = dPhi(i, kk) + (dNt(i, j) - Yt(i, j) * dLam[j]) * int_D(i, j) * exp(betaD * int_D(i, j)) -
                exp(betaD*int_D(i,j))*Yt(i,j)*dLam_dbetaD[j]-
                Yt(i, j-1) * int_dexpbetaD(i, j) * (D_status(i, j-1) * betaD + Covbeta[i]) -
                Yt(i, j-1) * int_expbetaD(i, j) * D_status(i, j-1);
            }
            else
            {
              dPhi(i, kk) = dPhi(i, kk) - Yt(i, j-1) * Covariates(i, kk-1) * int_expbetaD(i, j)-
                Yt(i, j) * exp(int_D(i, j)*betaD)*dLam_dbeta(kk-1, j);
            }
          }
        }
        
        
        
        
        
        // if (j == 1 && i ==0) arma::arma_print(res);
      }
      // printf("%d\n", i);
      for (int kk = 0; kk < p+1; kk++)
      {
        if (kk == 0)
        {
          fn[kk] = fn[kk] + IV_c[i] * res[i];
          
          for (int kkk = 0; kkk < p+1; kkk++)
          {
            Hessian(kk, kkk) = Hessian(kk, kkk) + IV_c[i] * dPhi(i, kkk);                        
          }
        }
        else
        {
          fn[kk] = fn[kk] + Covariates(i, kk-1) * res[i];
          for (int kkk = 0; kkk < p+1; kkk++)
          {
            // if(iter == 0) arma::arma_print(Covariates(i, kk-1) * dPhi(i, kkk));
            Hessian(kk, kkk) = Hessian(kk, kkk) + Covariates(i, kk-1) * dPhi(i, kkk);
          }
        }
      }
      
      
    }
    
    // arma::arma_print(Hessian);
    
    pk = -arma::solve(Hessian, fn);
    new_parameters = init_parameters + pk;
    fn_abs = arma::abs(fn);
    // arma::arma_print(fn_abs);
    // arma::arma_print(arma::sum(fn_abs));
    diff = arma::abs(new_parameters - init_parameters);
    if(arma::sum(fn_abs) < tol || arma::sum(diff) < (tol * tol))
    {
      step = iter;
      Convergence = true;
      break;
    }
    init_parameters = new_parameters;
    // prev_fn = fn;
    if (iter >= (max_iter - 1))
    {
      step = iter;
      Convergence = false;
    }
  }
  
  
  // for (int j = 0; j < k; j++)
  // {
  //   if (j == 0)
  //   {
  //     dLam[j] = 0.25*stime[j];
  //   } else 
  //   {
  //     dLam[j] = 0.25*(stime[j] - stime[j-1]);
  //   }
  // }
  
  
  
  
  // arma::arma_print(res);
  
  // if (iter > 0)
  // {
  //     pk_mat.col(0) = pk;
  //     H_com = pk_mat.t()*Hessian*pk_mat;
  //     H_com_val = H_com(0, 0);
  //     fn_2 = arma::pow(fn, 2);
  //     prev_fn_2 = arma::pow(prev_fn, 2);
  //     bktrkg = (arma::sum(fn_2)/2 - arma::sum(prev_fn_2)/2)/(arma::dot(pk, fn) + H_com_val);
  //     if(bktrkg > eta)
  //     {
  //         printf("here\n");
  //         init_parameters = init_parameters - (1 - contraction) * pk;
  //         if (iter >= max_iter - 1) Convergence = false;
  //         continue;
  //     }
  // }
  
  // arma::arma_print(fn);
  // arma::arma_print(Hessian);
  // arma::arma_print(dNt);
  // if(iter == 0)
  //     arma::arma_print(dPhi);
  // printf("This way\n");
  
  
  // asymptotic variance
  
  
  Asyvar = 0.0;
  hatsigma = 0.0;
  hatphi = 0.0;
  ave_hatsigma = 0.0;
  hatsigma_v.zeros();
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < k; j++)
    {
      if (j == 0)
      {
        ave_hatsigma = ave_hatsigma + (IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j) - Yt(i, j) * dLam[j]) - 
          IV_c[i] * int_expbetaD(i, j)* 1*(betaD* IV[i] + Covbeta[i]))/n;
      } else
      {
        ave_hatsigma = ave_hatsigma + (IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * dLam[j]) - 
          IV_c[i] * int_expbetaD(i, j)* Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]))/n;
      }
    }
  }
  
  
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < k; j++)
    {
      // hatsigma = hatsigma + ((IV_c[i] * IV_c[i])*exp(2*int_D(i, j)*betaD)*dNt(i, j));
      
      if(j==0)
      {
        hatsigma_v[i] = hatsigma_v[i] + (IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * dLam[j]) - 
          IV_c[i] * int_expbetaD(i, j)* 1*(betaD*IV[i] + Covbeta[i]));
        hatphi = hatphi + (IV_c[i] * int_D(i, j)* exp(int_D(i,j)*betaD)*(dNt(i, j) - Yt(i, j)*dLam[j]) -
          IV_c[i] * int_dexpbetaD(i, j)*1*(betaD*IV[i] + Covbeta[i]) -
          IV_c[i] * int_expbetaD(i, j)*1*IV[i]);
      } else 
      {
        hatsigma_v[i] = hatsigma_v[i] + IV_c[i]*exp(int_D(i, j) * betaD) *(dNt(i, j)- Yt(i, j) * dLam[j]) - 
          IV_c[i] * int_expbetaD(i, j)* Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]);
        hatphi = hatphi + (IV_c[i] * int_D(i, j)* exp(int_D(i,j)*betaD)*(dNt(i, j)- Yt(i, j)*dLam[j]) -
          IV_c[i] * int_dexpbetaD(i, j)*Yt(i, j-1)*(betaD*D_status(i, j-1) + Covbeta[i]) -
          IV_c[i] * int_expbetaD(i, j)*Yt(i, j-1)*D_status(i, j-1));
      }
      
    }
    hatsigma = hatsigma + (hatsigma_v[i] - ave_hatsigma)*(hatsigma_v[i] - ave_hatsigma);
  }
  Asyvar = hatsigma / (hatphi*hatphi);
  // arma::arma_print((hatsigma/ (Hessian(0, 0) * Hessian(0, 0))));
  
  
  return Rcpp::List::create(
    Named("x") = init_parameters,
    Named("Convergence") = Convergence,
    Named("iter") = step,
    Named("var") = Asyvar,
    Named("dLam") = dLam
  );
}
