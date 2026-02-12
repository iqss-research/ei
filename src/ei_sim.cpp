#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Forward declaration of bounds_cpp from bounds.cpp
NumericMatrix bounds_cpp(NumericVector x, NumericVector t, NumericVector n);

// Standard normal CDF
static double pnorm_fast(double x) {
  return 0.5 * erfc(-x * M_SQRT1_2);
}

// Bivariate normal CDF using Drezner & Wesolowsky (1990) algorithm
// P(X1 <= dh, X2 <= dk | correlation = r)
static double bvnorm_cdf(double dh, double dk, double r) {
  if (std::abs(r) < 1e-15) {
    return pnorm_fast(dh) * pnorm_fast(dk);
  }
  if (std::abs(r - 1.0) < 1e-15) {
    return pnorm_fast(std::min(dh, dk));
  }
  if (std::abs(r + 1.0) < 1e-15) {
    if (dh + dk < 0) return 0.0;
    return pnorm_fast(dh) + pnorm_fast(dk) - 1.0;
  }

  double tp = 2.0 * M_PI;
  double bvn = 0.0;

  if (std::abs(r) < 0.925) {
    static const double w6[] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910};
    static const double x6[] = {0.9324695142031521, 0.6612093864662645, 0.2386191860831969};

    double hs = (dh * dh + dk * dk) / 2.0;
    double asr = asin(r);

    for (int i = 0; i < 3; i++) {
      for (int is = -1; is <= 1; is += 2) {
        double sn = sin(asr * (is * x6[i] + 1.0) / 2.0);
        bvn += w6[i] * exp((sn * dh * dk - hs) / (1.0 - sn * sn));
      }
    }
    bvn *= asr / (2.0 * tp);
    bvn += pnorm_fast(-dh) * pnorm_fast(-dk);
  } else {
    if (r < 0) {
      dk = -dk;
    }
    double hk = dh * dk;
    double ass = (1.0 - r) * (1.0 + r);
    double a = sqrt(ass);
    double bs = (dh - dk) * (dh - dk);
    double c = (4.0 - hk) / 8.0;
    double d = (12.0 - hk) / 16.0;
    double asr = -(bs / ass + hk) / 2.0;

    if (asr > -100.0) {
      bvn = a * exp(asr) * (1.0 - c * (bs - ass) * (1.0 - d * bs / 5.0) / 3.0 + c * d * ass * ass / 5.0);
    }
    if (-hk < 100.0) {
      double b = sqrt(bs);
      bvn -= exp(-hk / 2.0) * sqrt(tp) * pnorm_fast(-b / a) * b *
        (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
    }
    a /= 2.0;

    static const double w10[] = {0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
                                  0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
                                  0.1316886384491766, 0.1420961093183820, 0.1491729864726037,
                                  0.1527533871307258};
    static const double x10[] = {0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
                                  0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
                                  0.5108670019508271, 0.3737060887154195, 0.2277858511416451,
                                  0.0765265211334973};

    for (int i = 0; i < 10; i++) {
      for (int is = -1; is <= 1; is += 2) {
        double xs = a * (is * x10[i] + 1.0);
        xs = xs * xs;
        double rs = sqrt(1.0 - xs);
        asr = -(bs / xs + hk) / 2.0;
        if (asr > -100.0) {
          bvn += a * w10[i] * exp(asr) * (exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
            (1.0 + c * xs * (1.0 + d * xs)));
        }
      }
    }
    bvn = -bvn / tp;

    if (r > 0) {
      bvn += pnorm_fast(-std::max(dh, dk));
    } else {
      bvn = -bvn;
      if (dk > dh) {
        bvn += pnorm_fast(dk) - pnorm_fast(dh);
      }
    }
  }
  return std::max(0.0, std::min(1.0, bvn));
}

// P(lower1 < X1 < upper1, lower2 < X2 < upper2 | correlation = rho)
static double pmvnorm_rect(double lower1, double upper1, double lower2, double upper2, double rho) {
  double p = bvnorm_cdf(upper1, upper2, rho);
  if (std::isfinite(lower1)) p -= bvnorm_cdf(lower1, upper2, rho);
  if (std::isfinite(lower2)) p -= bvnorm_cdf(upper1, lower2, rho);
  if (std::isfinite(lower1) && std::isfinite(lower2)) p += bvnorm_cdf(lower1, lower2, rho);
  return std::max(p, 1e-322);
}

// Vectorized createR for Rfun 1, 2, or 5
// [[Rcpp::export]]
NumericVector createR_cpp(NumericVector bb, NumericVector bw,
                          double sb, double sw, double rho,
                          LogicalVector sub, int Rfun) {
  int n = sub.length();
  std::vector<int> idx;
  for (int i = 0; i < n; i++) {
    if (sub[i]) idx.push_back(i);
  }
  int nsub = idx.size();
  NumericVector out(nsub);

  if (Rfun == 5) {
    double lower1 = -bb[idx[0]] / sb;
    double upper1 = lower1 + 1.0 / sb;
    double lower2 = -bw[idx[0]] / sw;
    double upper2 = lower2 + 1.0 / sw;
    double qi = pmvnorm_rect(lower1, upper1, lower2, upper2, rho);
    double logqi = log(qi);
    if (!R_finite(logqi) || R_IsNA(logqi)) logqi = 999.0;
    for (int i = 0; i < nsub; i++) out[i] = logqi;
    return out;
  }

  for (int i = 0; i < nsub; i++) {
    int j = idx[i];
    double lower1 = -bb[j] / sb;
    double upper1 = lower1 + 1.0 / sb;
    double lower2 = -bw[j] / sw;
    double upper2 = lower2 + 1.0 / sw;
    double qi = pmvnorm_rect(lower1, upper1, lower2, upper2, rho);
    double logqi = log(qi);
    if (!R_finite(logqi) || R_IsNA(logqi)) logqi = 999.0;
    out[i] = logqi;
  }
  return out;
}

// Full like() in C++ for the no-covariates case (Rfun==5, numb==1)
// [[Rcpp::export]]
double like_cpp(NumericVector param, NumericVector y, NumericVector x,
                NumericVector n_vec, int numb,
                double erho, double esigma, double ebeta) {

  double Bb0 = param[0];
  double Bw0 = param[1];
  double sb0 = param[2];
  double sw0 = param[3];
  double rho0 = param[4];

  double sb = exp(sb0);
  double sw = exp(sw0);
  double e2rho0 = exp(2.0 * rho0);
  double rho = (e2rho0 - 1.0) / (e2rho0 + 1.0);
  double sigb2 = sb * sb;
  double sigw2 = sw * sw;
  double sigbw = rho * sb * sw;

  double bb_val = Bb0 * (0.25 + sigb2) + 0.5;
  double bw_val = Bw0 * (0.25 + sigw2) + 0.5;

  int np = x.length();
  double enumtol = 0.0001;

  NumericMatrix bounds = bounds_cpp(x, y, n_vec);

  // createR for Rfun==5: compute once
  double lower1_R = -bb_val / sb;
  double upper1_R = lower1_R + 1.0 / sb;
  double lower2_R = -bw_val / sw;
  double upper2_R = lower2_R + 1.0 / sw;
  double R_val = log(pmvnorm_rect(lower1_R, upper1_R, lower2_R, upper2_R, rho));
  if (!R_finite(R_val) || R_IsNA(R_val)) R_val = 999.0;

  // Sigma for cT0/cT1 dmvnorm
  double det_sigma = sigb2 * sigw2 - sigbw * sigbw;
  double log_det_sigma = log(std::max(det_sigma, 1e-322));
  double log2pi2 = 2.0 * log(2.0 * M_PI);

  double llik_het = 0.0;
  double llik_wh = 0.0;
  double llik_bl = 0.0;
  double llik_cT0 = 0.0;
  double llik_cT1 = 0.0;

  for (int i = 0; i < np; i++) {
    double xi = x[i];
    double yi = y[i];

    if (xi == 0.0) {
      double epsilon = yi - bw_val;
      llik_wh += -0.5 * (log(sigw2) + (epsilon * epsilon) / sigw2);
      double Ebb = bb_val + rho * (sb / sw) * epsilon;
      double vbb = sigb2 * (1.0 - rho * rho);
      if (vbb < 1e-32) vbb = 0.0001;
      double s = sqrt(vbb);
      double b_s = (1.0 - Ebb) / s;
      double a_s = (0.0 - Ebb) / s;
      llik_wh += log(std::max(pnorm_fast(-a_s) - pnorm_fast(-b_s), 1e-322));
      llik_wh -= R_val;
    } else if (xi == 1.0) {
      double epsilon = yi - bb_val;
      llik_bl += -0.5 * (log(sigb2) + (epsilon * epsilon) / sigb2);
      double Ebb = bw_val + rho * (sw / sb) * epsilon;
      double vbb = sigw2 * (1.0 - rho * rho);
      if (vbb < 1e-32) vbb = 0.0001;
      double s = sqrt(vbb);
      double b_s = (1.0 - Ebb) / s;
      double a_s = (0.0 - Ebb) / s;
      llik_bl += log(std::max(pnorm_fast(-a_s) - pnorm_fast(-b_s), 1e-322));
      llik_bl -= R_val;
    } else if (yi < enumtol) {
      double d1 = 0.0 - bb_val;
      double d2 = 0.0 - bw_val;
      double maha = (sigw2 * d1 * d1 - 2.0 * sigbw * d1 * d2 + sigb2 * d2 * d2) / det_sigma;
      llik_cT0 += -0.5 * (log2pi2 + log_det_sigma + maha);
      llik_cT0 -= R_val;
    } else if (yi > 1.0 - enumtol) {
      double d1 = 1.0 - bb_val;
      double d2 = 1.0 - bw_val;
      double maha = (sigw2 * d1 * d1 - 2.0 * sigbw * d1 * d2 + sigb2 * d2 * d2) / det_sigma;
      llik_cT1 += -0.5 * (log2pi2 + log_det_sigma + maha);
      llik_cT1 -= R_val;
    } else {
      double omxi = 1.0 - xi;
      double s2 = sigb2 * xi * xi + sigw2 * omxi * omxi + 2.0 * sigbw * xi * omxi;
      double omega = sigb2 * xi + sigbw * omxi;
      double epsilon = yi - bb_val * xi - bw_val * omxi;
      double ebb = bb_val + (omega / s2) * epsilon;
      double vbb = sigb2 - (omega * omega) / s2;
      if (vbb < 1e-32) vbb = 0.0001;
      double s = sqrt(vbb);

      llik_het += -0.5 * (log(s2) + (epsilon * epsilon) / s2);
      double b_s = (bounds(i, 1) - ebb) / s;
      double a_s = (bounds(i, 0) - ebb) / s;
      llik_het += log(std::max(pnorm_fast(-a_s) - pnorm_fast(-b_s), 1e-322));
      llik_het -= R_val;
    }
  }

  double llik = llik_het + llik_bl + llik_wh + llik_cT0 + llik_cT1;

  // Priors
  double prior = 0.0;
  double lpdfnorm = -0.5 * log(2.0 * M_PI) - log(erho) - 0.5 * (rho0 * rho0) / (erho * erho);
  if (esigma > 0) prior -= (1.0 / (2.0 * esigma * esigma)) * (sigb2 + sigw2);
  if (erho > 0) prior += lpdfnorm;
  if (ebeta > 0 && bb_val < 0) prior -= 0.5 * (bb_val * bb_val / ebeta);
  if (ebeta > 0 && bb_val > 1) prior -= 0.5 * ((bb_val - 1.0) * (bb_val - 1.0) / ebeta);
  if (ebeta > 0 && bw_val < 0) prior -= 0.5 * (bw_val * bw_val / ebeta);
  if (ebeta > 0 && bw_val > 1) prior -= 0.5 * ((bw_val - 1.0) * (bw_val - 1.0) / ebeta);

  llik += prior;
  if (R_IsNA(llik) || !R_finite(llik)) llik = -1e20;
  return -llik;
}

// Batch computation of like_cpp for importance sampling
// [[Rcpp::export]]
NumericVector like_batch_cpp(NumericMatrix draw, NumericVector y, NumericVector x,
                             NumericVector n_vec, int numb,
                             double erho, double esigma, double ebeta) {
  int nsims = draw.nrow();
  NumericVector result(nsims);
  for (int i = 0; i < nsims; i++) {
    NumericVector param = draw(i, _);
    result[i] = -like_cpp(param, y, x, n_vec, numb, erho, esigma, ebeta);
  }
  return result;
}
