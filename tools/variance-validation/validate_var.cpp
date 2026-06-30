// Standalone validation of the DRIV.s variance estimators.
//
// Compares, by Monte Carlo:
//   - empirical SD of theta-hat across replicates  (the truth)
//   - mean of the ORIGINAL scalar-sandwich SE      (hatsigma / hatphi^2)
//   - mean of the NEW joint-sandwich SE            (.tex semi-parametric section)
// plus 95% CI coverage for each SE.
//
// The estimator core is a faithful port of src/driv_s_est.cpp (Rcpp stripped),
// with the j==0 dPhi indexing bug IV[j] -> IV[i] fixed.
//
// build:  g++ -O2 validate_var.cpp -o validate_var -larmadillo
#include <armadillo>
#include <cstdio>
#include <random>
using namespace arma;

static inline double expit(double d){ return 1.0/(1.0+std::exp(-d)); }

struct FitOut {
  vec    x;          // (theta, alpha)
  bool   conv;
  double var_orig;   // original scalar-sandwich Var(theta_hat)
  double var_new;    // new joint-sandwich  Var(theta_hat)
};

// Zdesign: n x pp propensity design (incl. intercept) used to form IV_c and the
// gamma block.  IV_c[i] = IV[i] - expit(Zdesign_i . gamma_hat) is precomputed
// outside and passed in; Zdesign + IV_c let us reconstruct expit_i = IV-IV_c.
FitOut fit(vec init_parameters, const vec& time, const vec& event, const vec& IV,
           const vec& IV_c, const mat& Covariates, const mat& D_status,
           const vec& stime, const mat& Zdesign,
           int max_iter, double tol)
{
  int n = time.n_elem;
  int k = stime.n_elem;
  int p = init_parameters.n_elem - 1;
  int pp = Zdesign.n_cols;          // propensity dim (incl intercept)

  double betaD=0;
  bool Convergence=false; int /*step=0,*/ used_iter=0;
  vec beta(p), res(n), Covbeta(n), fn(p+1), SY(k), dSY(k), dLam_dbetaD(k), dLam(k), pk(p+1), new_parameters(p+1);
  mat int_D(n,k), dNt(n,k), Yt(n,k), int_expbetaD(n,k), int_dexpbetaD(n,k), dLam_dbeta(p,k), dPhi(n,p+1), Hessian(p+1,p+1);

  for (int iter=0; iter<max_iter; iter++)
  {
    res.zeros(); Covbeta.zeros(); dPhi.zeros(); Hessian.zeros();
    fn.zeros(); SY.zeros(); dSY.zeros(); int_D.zeros(); int_expbetaD.zeros(); int_dexpbetaD.zeros();
    betaD = init_parameters[0];
    beta  = init_parameters.subvec(1,p);

    for (int i=0;i<n;i++) for (int ii=0;ii<p;ii++) Covbeta[i]+=Covariates(i,ii)*beta[ii];

    for (int i=0;i<n;i++) for (int j=0;j<k;j++){
      if(j==0) int_D(i,j)=IV[i]*stime[j];
      else     int_D(i,j)=int_D(i,j-1)+D_status(i,j-1)*(stime[j]-stime[j-1]);
      dNt(i,j)=((time[i]==stime[j])?1:0)*event[i];
      Yt(i,j)=(time[i]>=stime[j])?1:0;
      SY[j]+=std::exp(betaD*int_D(i,j))*Yt(i,j);
      dSY[j]+=int_D(i,j)*std::exp(betaD*int_D(i,j))*Yt(i,j);
    }

    dLam.zeros(); dLam_dbetaD.zeros(); dLam_dbeta.zeros();
    for (int i=0;i<n;i++) for (int j=0;j<k;j++){
      double dt = (j==0)? stime[j] : (stime[j]-stime[j-1]);
      double e = std::exp(int_D(i,j)*betaD);
      dLam[j]      += Yt(i,j)*e*(dNt(i,j)-D_status(i,j)*betaD*dt - Covbeta[i]*dt)/SY[j];
      dLam_dbetaD[j]+= Yt(i,j)*int_D(i,j)*e*(dNt(i,j)-(D_status(i,j)*betaD+Covbeta[i])*dt)/SY[j]
                     - Yt(i,j)*e*D_status(i,j)*dt/SY[j]
                     - dSY[j]*Yt(i,j)*e*(dNt(i,j)-D_status(i,j)*betaD*dt-Covbeta[i]*dt)/(SY[j]*SY[j]);
      for (int kk=0;kk<p;kk++)
        dLam_dbeta(kk,j) += -Yt(i,j)*e*Covariates(i,kk)*dt/SY[j];
    }

    for (int i=0;i<n;i++){
      for (int j=0;j<k;j++){
        double e = std::exp(betaD*int_D(i,j));
        if(j==0){
          if(betaD!=0){
            int_expbetaD(i,j)=(IV[i]>0.5)?(std::exp(betaD*int_D(i,j))-1)/betaD:stime[j];
            int_dexpbetaD(i,j)=(IV[i]>0.5)?(int_D(i,j)*std::exp(betaD*int_D(i,j))/betaD-int_expbetaD(i,j)/betaD):0;
          } else { int_expbetaD(i,j)=stime[j]; int_dexpbetaD(i,j)=int_D(i,j)*int_D(i,j)/2; }
          res[i]+= e*(dNt(i,j)-Yt(i,j)*dLam[j]) - int_expbetaD(i,j)*(IV[i]*betaD+Covbeta[i]);
          for(int kk=0;kk<p+1;kk++){
            if(kk==0)
              dPhi(i,kk)+= (dNt(i,j)-Yt(i,j)*dLam[j])*int_D(i,j)*e
                          - e*Yt(i,j)*dLam_dbetaD[j]
                          - int_dexpbetaD(i,j)*(IV[i]*betaD+Covbeta[i])   // fixed: IV[i]
                          - int_expbetaD(i,j)*IV[i];                       // fixed: IV[i]
            else
              dPhi(i,kk)+= -Covariates(i,kk-1)*int_expbetaD(i,j)
                          - Yt(i,j)*e*dLam_dbeta(kk-1,j);
          }
        } else {
          if(betaD!=0){
            int_expbetaD(i,j)=(D_status(i,j-1)>0.5)?(std::exp(betaD*int_D(i,j))-std::exp(betaD*int_D(i,j-1)))/betaD:(stime[j]-stime[j-1])*std::exp(betaD*int_D(i,j));
            int_dexpbetaD(i,j)=(D_status(i,j-1)>0.5)?((int_D(i,j)*std::exp(betaD*int_D(i,j))-int_D(i,j-1)*std::exp(betaD*int_D(i,j-1)))/betaD-int_expbetaD(i,j)/betaD):(int_D(i,j)*std::exp(betaD*int_D(i,j))*(stime[j]-stime[j-1]));
          } else {
            int_expbetaD(i,j)=stime[j]-stime[j-1];
            int_dexpbetaD(i,j)=(stime[j]-stime[j-1])*int_D(i,j-1)+D_status(i,j-1)*(stime[j]-stime[j-1])*(stime[j]-stime[j-1])/2;
          }
          res[i]+= e*(dNt(i,j)-Yt(i,j)*dLam[j]) - Yt(i,j-1)*int_expbetaD(i,j)*(D_status(i,j-1)*betaD+Covbeta[i]);
          for(int kk=0;kk<p+1;kk++){
            if(kk==0)
              dPhi(i,kk)+= (dNt(i,j)-Yt(i,j)*dLam[j])*int_D(i,j)*e
                          - e*Yt(i,j)*dLam_dbetaD[j]
                          - Yt(i,j-1)*int_dexpbetaD(i,j)*(D_status(i,j-1)*betaD+Covbeta[i])
                          - Yt(i,j-1)*int_expbetaD(i,j)*D_status(i,j-1);
            else
              dPhi(i,kk)+= -Yt(i,j-1)*Covariates(i,kk-1)*int_expbetaD(i,j)
                          - Yt(i,j)*e*dLam_dbeta(kk-1,j);
          }
        }
      }
      for(int kk=0;kk<p+1;kk++){
        double w = (kk==0)? IV_c[i] : Covariates(i,kk-1);
        fn[kk]+= w*res[i];
        for(int kkk=0;kkk<p+1;kkk++) Hessian(kk,kkk)+= w*dPhi(i,kkk);
      }
    }

    pk = -arma::solve(Hessian, fn);
    new_parameters = init_parameters + pk;
    used_iter=iter;
    if(arma::sum(arma::abs(fn))<tol || arma::sum(arma::abs(new_parameters-init_parameters))<(tol*tol)){
      Convergence=true; init_parameters=new_parameters; break;
    }
    init_parameters=new_parameters;
    if(iter>=max_iter-1) Convergence=false;
  }
  betaD=init_parameters[0]; beta=init_parameters.subvec(1,p);

  // ---------- ORIGINAL scalar-sandwich variance ----------
  double ave_hatsigma=0, hatsigma=0, hatphi=0;
  vec hatsigma_v(n, fill::zeros), res_store(n, fill::zeros);
  for(int i=0;i<n;i++) for(int j=0;j<k;j++){
    double e=std::exp(int_D(i,j)*betaD);
    double term;
    if(j==0) term = IV_c[i]*e*(dNt(i,j)-Yt(i,j)*dLam[j]) - IV_c[i]*int_expbetaD(i,j)*1*(betaD*IV[i]+Covbeta[i]);
    else     term = IV_c[i]*e*(dNt(i,j)-Yt(i,j)*dLam[j]) - IV_c[i]*int_expbetaD(i,j)*Yt(i,j-1)*(betaD*D_status(i,j-1)+Covbeta[i]);
    ave_hatsigma += term/n;
  }
  for(int i=0;i<n;i++){
    for(int j=0;j<k;j++){
      double e=std::exp(int_D(i,j)*betaD);
      if(j==0){
        hatsigma_v[i]+= IV_c[i]*e*(dNt(i,j)-Yt(i,j)*dLam[j]) - IV_c[i]*int_expbetaD(i,j)*1*(betaD*IV[i]+Covbeta[i]);
        hatphi += IV_c[i]*int_D(i,j)*e*(dNt(i,j)-Yt(i,j)*dLam[j])
                - IV_c[i]*int_dexpbetaD(i,j)*1*(betaD*IV[i]+Covbeta[i])
                - IV_c[i]*int_expbetaD(i,j)*1*IV[i];
      } else {
        hatsigma_v[i]+= IV_c[i]*e*(dNt(i,j)-Yt(i,j)*dLam[j]) - IV_c[i]*int_expbetaD(i,j)*Yt(i,j-1)*(betaD*D_status(i,j-1)+Covbeta[i]);
        hatphi += IV_c[i]*int_D(i,j)*e*(dNt(i,j)-Yt(i,j)*dLam[j])
                - IV_c[i]*int_dexpbetaD(i,j)*Yt(i,j-1)*(betaD*D_status(i,j-1)+Covbeta[i])
                - IV_c[i]*int_expbetaD(i,j)*Yt(i,j-1)*D_status(i,j-1);
      }
    }
    hatsigma += (hatsigma_v[i]-ave_hatsigma)*(hatsigma_v[i]-ave_hatsigma);
  }
  double var_orig = hatsigma/(hatphi*hatphi);

  // ---------- NEW joint-sandwich variance ----------
  // risk-set weighted averages xi_IV[j], xi_X(m,j)
  vec SY_IVc(k, fill::zeros); mat SY_X(p,k, fill::zeros);
  for(int i=0;i<n;i++) for(int j=0;j<k;j++){
    double w=std::exp(betaD*int_D(i,j))*Yt(i,j);
    SY_IVc[j]+= w*IV_c[i];
    for(int m=0;m<p;m++) SY_X(m,j)+= w*Covariates(i,m);
  }
  int d = 1 + p + pp;
  mat U(d,n, fill::zeros);              // Utilde per subject (columns)
  for(int i=0;i<n;i++){
    double res_i=0, BPsi=0; vec BPhi(p, fill::zeros);
    for(int j=0;j<k;j++){
      double e=std::exp(int_D(i,j)*betaD);
      double dres;
      if(j==0) dres = e*(dNt(i,j)-Yt(i,j)*dLam[j]) - int_expbetaD(i,j)*1*(betaD*IV[i]+Covbeta[i]);
      else     dres = e*(dNt(i,j)-Yt(i,j)*dLam[j]) - Yt(i,j-1)*int_expbetaD(i,j)*(betaD*D_status(i,j-1)+Covbeta[i]);
      res_i += dres;
      double xiIV = (SY[j]>0)? SY_IVc[j]/SY[j] : 0;
      BPsi += xiIV*dres;
      for(int m=0;m<p;m++){ double xiX=(SY[j]>0)?SY_X(m,j)/SY[j]:0; BPhi[m]+=xiX*dres; }
    }
    res_store[i]=res_i;
    BPsi=-BPsi; for(int m=0;m<p;m++) BPhi[m]=-BPhi[m];
    U(0,i)= IV_c[i]*res_i + BPsi;                          // Psi + B_Psi
    for(int m=0;m<p;m++)  U(1+m,i)= Covariates(i,m)*res_i + BPhi[m];   // Phi + B_Phi
    for(int m=0;m<pp;m++) U(1+p+m,i)= Zdesign(i,m)*IV_c[i];           // propensity block (no corr)
  }
  vec Ubar = mean(U,1);
  for(int i=0;i<n;i++) U.col(i)-=Ubar;
  mat Omega = U*U.t();                  // sum of centered outer products

  // bread A (sums)
  mat A(d,d, fill::zeros);
  A.submat(0,0,p,p) = Hessian;          // d(Psi,Phi)/d(theta,alpha), profiled
  for(int i=0;i<n;i++){
    double ex=IV[i]-IV_c[i]; double w=ex*(1-ex);
    for(int m=0;m<pp;m++) A(0,1+p+m) += -w*res_store[i]*Zdesign(i,m);  // dPsi/dgamma
    for(int m1=0;m1<pp;m1++) for(int m2=0;m2<pp;m2++)
      A(1+p+m1,1+p+m2) += -w*Zdesign(i,m1)*Zdesign(i,m2);             // dg/dgamma
  }
  mat Ainv = arma::inv(A);
  mat V = Ainv*Omega*Ainv.t();
  double var_new = V(0,0);

  // ===== verbatim copy of the package block (src/driv_s_est.cpp) -> var_pkg =====
  double var_pkg;
  {
    int pp = p + 1;
    int d  = 1 + p + pp;
    arma::vec SY_IVc(k, arma::fill::zeros);
    arma::mat SY_X(p, k, arma::fill::zeros);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < k; j++) {
        double w = exp(betaD * int_D(i, j)) * Yt(i, j);
        SY_IVc[j] += w * IV_c[i];
        for (int m = 0; m < p; m++) SY_X(m, j) += w * Covariates(i, m);
      }
    arma::mat Ujoint(d, n, arma::fill::zeros);
    arma::vec res_store2(n, arma::fill::zeros);
    for (int i = 0; i < n; i++) {
      double res_i = 0.0, BPsi = 0.0;
      arma::vec BPhi(p, arma::fill::zeros);
      for (int j = 0; j < k; j++) {
        double e = exp(int_D(i, j) * betaD);
        double dres;
        if (j == 0)
          dres = e * (dNt(i, j) - Yt(i, j) * dLam[j]) - int_expbetaD(i, j) * (betaD * IV[i] + Covbeta[i]);
        else
          dres = e * (dNt(i, j) - Yt(i, j) * dLam[j]) - Yt(i, j-1) * int_expbetaD(i, j) * (betaD * D_status(i, j-1) + Covbeta[i]);
        res_i += dres;
        double xiIV = (SY[j] > 0) ? SY_IVc[j] / SY[j] : 0.0;
        BPsi += xiIV * dres;
        for (int m = 0; m < p; m++) { double xiX = (SY[j] > 0) ? SY_X(m, j) / SY[j] : 0.0; BPhi[m] += xiX * dres; }
      }
      res_store2[i] = res_i;
      Ujoint(0, i) = IV_c[i] * res_i - BPsi;
      for (int m = 0; m < p; m++) Ujoint(1 + m, i) = Covariates(i, m) * res_i - BPhi[m];
      Ujoint(1 + p, i) = IV_c[i];
      for (int m = 0; m < p; m++) Ujoint(1 + p + 1 + m, i) = Covariates(i, m) * IV_c[i];
    }
    arma::vec Ubar = arma::mean(Ujoint, 1);
    for (int i = 0; i < n; i++) Ujoint.col(i) -= Ubar;
    arma::mat Omega = Ujoint * Ujoint.t();
    arma::mat Aj(d, d, arma::fill::zeros);
    Aj.submat(0, 0, p, p) = Hessian;
    for (int i = 0; i < n; i++) {
      double ex = IV[i] - IV_c[i];
      double w  = ex * (1.0 - ex);
      Aj(0, 1 + p) += -w * res_store2[i];
      for (int m = 0; m < p; m++) Aj(0, 1 + p + 1 + m) += -w * res_store2[i] * Covariates(i, m);
      for (int a = 0; a < pp; a++) {
        double za = (a == 0) ? 1.0 : Covariates(i, a-1);
        for (int b = 0; b < pp; b++) {
          double zb = (b == 0) ? 1.0 : Covariates(i, b-1);
          Aj(1 + p + a, 1 + p + b) += -w * za * zb;
        }
      }
    }
    arma::mat Ainv;
    if (arma::inv(Ainv, Aj)) { arma::mat V = Ainv * Omega * Ainv.t(); var_pkg = V(0, 0); }
    else { var_pkg = var_orig; }
  }
  if (std::fabs(var_pkg - var_new) > 1e-9*std::fabs(var_new)+1e-12) {
    fprintf(stderr, "TRANSCRIPTION MISMATCH: var_new=%.8e var_pkg=%.8e\n", var_new, var_pkg);
  }

  FitOut out; out.x=init_parameters; out.conv=Convergence;
  out.var_orig=var_orig; out.var_new=var_new;
  return out;
}

// logistic regression (Newton) of y on design Z (n x pp), returns coef
vec logit_fit(const mat& Z, const vec& y, int iters=50){
  int pp=Z.n_cols; vec g(pp, fill::zeros);
  for(int it=0; it<iters; it++){
    vec eta=Z*g; vec mu(y.n_elem);
    for(uword i=0;i<mu.n_elem;i++) mu[i]=expit(eta[i]);
    vec wvec = mu%(1-mu);
    mat ZtWZ = Z.t()*(Z.each_col()%wvec);
    vec score = Z.t()*(y-mu);
    vec step;
    if(!arma::solve(step, ZtWZ, score)) break;
    g+=step;
    if(arma::norm(step,2)<1e-9) break;
  }
  return g;
}

int main(int argc, char** argv){
  int    R    = (argc>1)? std::atoi(argv[1]) : 500;   // replicates
  int    n    = (argc>2)? std::atoi(argv[2]) : 800;   // sample size
  int    scen = (argc>3)? std::atoi(argv[3]) : 0;     // 0=no switch, 1=switching
  int    p    = (argc>4)? std::atoi(argv[4]) : 1;     // number of covariates
  double theta= 0.10, lam0=0.25;
  vec alpha(p), gamma(p);
  for(int m=0;m<p;m++){ alpha[m]=(scen==1)?0.40:0.25; gamma[m]=(scen==1)?1.5:1.0; }
  double diffcoef=0.5, beta_sw=0.25;   // switching params
  double Tmax = 5.0, cens=0.06, gridstep=0.1;
  const double INF = std::numeric_limits<double>::infinity();

  std::vector<double> ths, se_o, se_n; ths.reserve(R); se_o.reserve(R); se_n.reserve(R);
  int conv_ct=0, cov_o=0, cov_n=0, valid=0;

  for(int r=0;r<R;r++){
    std::mt19937_64 gen(1000+r);
    std::uniform_real_distribution<double> U01(0,1);
    std::exponential_distribution<double> E1(1.0);

    mat X(n,p); vec Z(n), W(n), time(n), event(n);
    for(int i=0;i<n;i++) for(int m=0;m<p;m++) X(i,m)=U01(gen);
    // centered linear propensity
    vec Xg(n); double mxg=0;
    for(int i=0;i<n;i++){ double s=0; for(int m=0;m<p;m++) s+=gamma[m]*X(i,m); Xg[i]=s; mxg+=s/n; }
    for(int i=0;i<n;i++){ double pz=expit(Xg[i]-mxg); Z[i]=(U01(gen)<pz)?1.0:0.0; }
    for(int i=0;i<n;i++){
      double Xa=0, Xb=0; for(int m=0;m<p;m++){ Xa+=alpha[m]*X(i,m); Xb+=beta_sw*X(i,m); }
      // switching time W (scen==1), else never
      double Wi = INF;
      if(scen==1){
        double Wexp=E1(gen);
        double denom = 0.5*(0.1 + (Z[i]>0.5? diffcoef : -diffcoef) + Xb);
        Wi = (denom>0)? Wexp/denom : INF;
        if(Wi<=0) Wi=INF;
      }
      // additive-hazard survival with switching adjustment (Settings.R SurvTime)
      double Te=E1(gen);
      double lamD = lam0 + theta*Z[i] + Xa;
      double T_D = Te/lamD;
      if(std::isfinite(Wi) && T_D>=Wi){
        double lamSw = lam0 + theta*(1-Z[i]) + Xa;
        T_D = (T_D*lamD - (theta*Z[i]*Wi - theta*(1-Z[i])*Wi))/lamSw;
      }
      W[i]=Wi;
      double Ci=E1(gen)/cens;
      double td=std::min(T_D,Ci); td=std::min(td,Tmax);
      time[i]=std::ceil(td/gridstep)*gridstep;
      event[i]=(T_D<=Ci && T_D<=Tmax)?1.0:0.0;
    }
    // grid = sorted unique observed times
    vec ut = arma::unique(time); vec stime = arma::sort(ut);
    int k=stime.n_elem;
    // D_status(i,j): Z before switch, 1-Z after (no switch if W=Inf)
    mat D_status(n,k);
    for(int i=0;i<n;i++) for(int j=0;j<k;j++)
      D_status(i,j) = (stime[j] <= W[i]) ? Z[i] : (1.0-Z[i]);

    // propensity design with intercept; gamma_hat via logistic regression
    mat Zd(n,p+1); for(int i=0;i<n;i++){ Zd(i,0)=1.0; for(int m=0;m<p;m++) Zd(i,1+m)=X(i,m); }
    vec gh = logit_fit(Zd, Z);
    vec IV_c(n); for(int i=0;i<n;i++) IV_c[i]=Z[i]-expit(dot(Zd.row(i),gh.t()));

    vec init(p+1, fill::zeros);
    FitOut f = fit(init, time, event, Z, IV_c, X, D_status, stime, Zd, 50, 1e-6);
    if(!f.conv) continue;
    if(f.var_orig<=0 || f.var_new<=0 || !std::isfinite(f.var_orig) || !std::isfinite(f.var_new)) continue;
    conv_ct++; valid++;
    double th=f.x[0]; ths.push_back(th);
    double so=std::sqrt(f.var_orig), sn=std::sqrt(f.var_new);
    se_o.push_back(so); se_n.push_back(sn);
    if(std::fabs(th-theta)<=1.96*so) cov_o++;
    if(std::fabs(th-theta)<=1.96*sn) cov_n++;
  }

  // summaries
  auto mean=[&](std::vector<double>&v){ double s=0; for(double x:v)s+=x; return s/v.size(); };
  double mth=mean(ths);
  double sd_emp=0; for(double x:ths) sd_emp+=(x-mth)*(x-mth); sd_emp=std::sqrt(sd_emp/(ths.size()-1));
  double mse_o=mean(se_o), mse_n=mean(se_n);

  printf("=== DRIV.s variance validation (no-switching, additive hazard) ===\n");
  printf("R(valid)=%d  n=%d  theta_true=%.3f\n", valid, n, theta);
  printf("mean theta_hat        : %.4f   (bias %.4f)\n", mth, mth-theta);
  printf("empirical SD(theta)   : %.4f   <-- truth\n", sd_emp);
  printf("mean SE original      : %.4f   (ratio to truth %.3f)\n", mse_o, mse_o/sd_emp);
  printf("mean SE new (joint)   : %.4f   (ratio to truth %.3f)\n", mse_n, mse_n/sd_emp);
  printf("95%% CI coverage orig  : %.3f\n", (double)cov_o/valid);
  printf("95%% CI coverage new   : %.3f\n", (double)cov_n/valid);
  // per-replicate divergence between the two SE estimators
  double mad=0,mx=0; for(size_t i=0;i<se_o.size();i++){ double rd=std::fabs(se_n[i]-se_o[i])/se_o[i]; mad+=rd/se_o.size(); if(rd>mx)mx=rd; }
  printf("per-rep |SE_new-SE_orig|/SE_orig: mean %.4f  max %.4f\n", mad, mx);
  printf("first reps (theta, SE_orig, SE_new):\n");
  for(int i=0;i<5 && i<(int)ths.size();i++) printf("   %.4f  %.5f  %.5f\n", ths[i], se_o[i], se_n[i]);
  return 0;
}
