#pragma once

// =======================================================================
//  sqg_tl.hpp  -  tangent linear and adjoint of the SQG model (sqg.hpp)
//
//  SQGTL linearizes one "model step" of sqg.hpp as ijedi::ModelSQG runs it:
//
//      grid PV x  ->  q = rfft2(x)  ->  nsteps RK4 steps  ->  irfft2  ->  grid PV
//
//  about a trajectory given by the grid PV at the start of the step. The trajectory
//  through the RK4 stages is recomputed in double precision from that grid; TL and
//  adjoint run in double precision too, so the adjoint is exact to rounding (the
//  nonlinear sqg.hpp model runs in single precision).
//
//  The RK4 scheme is sqg.hpp's: the hyperdiffusion term uses the state at the start
//  of the step in all four stages. The only nonlinearity is the Jacobian (advection)
//  term; everything else -- FFTs, PV inversion, zero-padding / truncation for the
//  dealiased Jacobian, relaxation, Ekman damping, hyperdiffusion -- is (real-)linear.
//
//  Adjoints are taken with respect to the Euclidean inner product of the grid values,
//  i.e. that of JEDI increments (oops::dotProductFieldSets). Every complex-to-real
//  inverse FFT first projects its input onto Hermitian symmetry (in the k=0 and Nyquist
//  columns), so the TL/AD pair is consistent whatever an FFT library does with
//  non-Hermitian input; with MKL's FFTW interface, c2r already acts that way, so the TL
//  is the exact derivative of the nonlinear model.
//
//  Dependencies: FFTW3 double precision (fftw_*) or Intel MKL's FFTW3 interface.
// =======================================================================

#include <fftw3.h>

#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

class SQGTL {
 public:
  using cplx = std::complex<double>;
  using RVec = std::vector<double>;
  using CVec = std::vector<cplx>;

  SQGTL(int N, double f, double nsq, double L, double H, double U, double r, double tdiab,
        int diff_order, double diff_efold, double dt)
      : N_(N), Nc_(N / 2 + 1), Np_(3 * N / 2), Ncp_(3 * N / 4 + 1),
        tdiab_(tdiab), dt_(dt), fft_(N), fftp_(3 * N / 2) {
    if (N % 2) throw std::invalid_argument("SQGTL: N must be even");
    r_[0] = r; r_[1] = -r;
    const double pi = M_PI;
    const int M = N_ * Nc_, Mp = Np_ * Ncp_;
    ksq_.resize(M); hovermu_.resize(M); tanhmu_.resize(M); sinhmu_.resize(M); hyperdiff_.resize(M);
    const double ktotcutoff = pi * N_ / L;
    const double eps = std::numeric_limits<float>::epsilon();  // as sqg.hpp
    for (int i = 0; i < N_; ++i) {
      const double lv = 2. * pi * ((i <= N_ / 2) ? i : i - N_) / L;
      for (int j = 0; j < Nc_; ++j) {
        const int idx = i * Nc_ + j;
        const double kv = 2. * pi * j / L;
        ksq_[idx] = kv * kv + lv * lv;
        double mu = std::sqrt(ksq_[idx]) * std::sqrt(nsq) * H / f;
        if (mu < eps) mu = eps;
        hovermu_[idx] = H / mu;
        tanhmu_[idx] = std::tanh(mu);
        sinhmu_[idx] = std::sinh(mu);
        hyperdiff_[idx] = -(1. / diff_efold) * std::pow(std::sqrt(ksq_[idx]) / ktotcutoff,
                                                        double(diff_order));
      }
    }
    ikp_.resize(Mp); ilp_.resize(Mp);
    for (int i = 0; i < Np_; ++i) {
      const double lv = 2. * pi * ((i <= Np_ / 2) ? i : i - Np_) / L;
      for (int j = 0; j < Ncp_; ++j) {
        ikp_[i * Ncp_ + j] = cplx(0., 2. * pi * j / L);
        ilp_[i * Ncp_ + j] = cplx(0., lv);
      }
    }
    // basic-state PV (relaxation target), as sqg.hpp
    const double l = 2. * pi / L, mu0 = l * std::sqrt(nsq) * H / f;
    RVec pvbar(2 * N_ * N_);
    for (int k = 0; k < 2; ++k)
      for (int j = 0; j < N_; ++j)
        for (int i = 0; i < N_; ++i)
          pvbar[k * N_ * N_ + j * N_ + i] = -(mu0 * 0.5 * U / (l * H)) * std::cosh(0.5 * mu0)
                                            * std::cos(l * j * L / N_) / std::sinh(0.5 * mu0);
    pveq_ = R(pvbar);
  }

  /// Nonlinear model step in double precision (for testing): grid -> grid.
  RVec nl(const RVec & x, int nsteps) const {
    CVec q = R(x);
    for (int n = 0; n < nsteps; ++n) q = rk4(q, nullptr);
    return scaleC(q);
  }

  /// Tangent linear of one model step about trajectory grid xtraj: dx -> dy.
  RVec tl(const RVec & xtraj, const RVec & dx, int nsteps) const {
    CVec q = R(xtraj), dq = R(dx);
    CVec stages[4];
    for (int n = 0; n < nsteps; ++n) {
      CVec qn = rk4(q, stages);
      const double a = 0.5 * dt_;
      CVec dk1 = tltend(dq, dq, stages[0]);
      CVec dk2 = tltend(axpy(dq, a, dk1), dq, stages[1]);
      CVec dk3 = tltend(axpy(dq, a, dk2), dq, stages[2]);
      CVec dk4 = tltend(axpy(dq, dt_, dk3), dq, stages[3]);
      for (size_t i = 0; i < dq.size(); ++i)
        dq[i] += dt_ / 6. * (dk1[i] + 2. * dk2[i] + 2. * dk3[i] + dk4[i]);
      q = qn;
    }
    return scaleC(dq);
  }

  /// Adjoint of tl(): dy -> dx (same trajectory).
  RVec ad(const RVec & xtraj, const RVec & dy, int nsteps) const {
    // forward sweep: trajectory states of all RK4 stages
    std::vector<CVec> traj(4 * nsteps);
    CVec q = R(xtraj);
    for (int n = 0; n < nsteps; ++n) q = rk4(q, &traj[4 * n]);
    // adjoint of dy = normN * C(dq)
    const double normN = 1. / (double(N_) * N_);
    RVec ys(dy);
    for (double & v : ys) v *= normN;
    CVec dqa = CT(ys);
    const double a = 0.5 * dt_;
    for (int n = nsteps - 1; n >= 0; --n) {
      const CVec * st = &traj[4 * n];
      CVec dk1a = scale(dqa, dt_ / 6.), dk2a = scale(dqa, dt_ / 3.);
      CVec dk3a = scale(dqa, dt_ / 3.), dk4a = scale(dqa, dt_ / 6.);
      CVec qna(dqa);  // the "dq_n +" term
      CVec sa, na;
      adtend(dk4a, st[3], sa, na);            // stage 4: dq4 = dq_n + dt*dk3
      add(qna, na); add(qna, sa); add(dk3a, scale(sa, dt_));
      adtend(dk3a, st[2], sa, na);            // stage 3: dq3 = dq_n + a*dk2
      add(qna, na); add(qna, sa); add(dk2a, scale(sa, a));
      adtend(dk2a, st[1], sa, na);            // stage 2: dq2 = dq_n + a*dk1
      add(qna, na); add(qna, sa); add(dk1a, scale(sa, a));
      adtend(dk1a, st[0], sa, na);            // stage 1: state dq_n
      add(qna, na); add(qna, sa);
      dqa = qna;
    }
    return RT(dqa);
  }

 private:
  // --------------------------------------------------------------- FFTs (one level)
  struct Fft {
    int n, nc;
    double * rbuf;
    fftw_complex * cbuf;
    fftw_plan fwd, bwd;
    explicit Fft(int n_) : n(n_), nc(n_ / 2 + 1) {
      rbuf = static_cast<double *>(fftw_malloc(sizeof(double) * n * n));
      cbuf = static_cast<fftw_complex *>(fftw_malloc(sizeof(fftw_complex) * n * nc));
      if (!rbuf || !cbuf) throw std::bad_alloc();
      fwd = fftw_plan_dft_r2c_2d(n, n, rbuf, cbuf, FFTW_ESTIMATE);
      bwd = fftw_plan_dft_c2r_2d(n, n, cbuf, rbuf, FFTW_ESTIMATE);
    }
    ~Fft() { fftw_destroy_plan(fwd); fftw_destroy_plan(bwd); fftw_free(rbuf); fftw_free(cbuf); }
    Fft(const Fft &) = delete;
    Fft & operator=(const Fft &) = delete;
    void r2c(const double * in, cplx * out) const {
      std::memcpy(rbuf, in, sizeof(double) * n * n);
      fftw_execute_dft_r2c(fwd, rbuf, cbuf);
      std::memcpy(out, cbuf, sizeof(cplx) * n * nc);
    }
    // unnormalized inverse of a Hermitian-projected input (k=0 and Nyquist columns)
    void c2r(const cplx * in, double * out) const {
      cplx * c = reinterpret_cast<cplx *>(cbuf);
      std::memcpy(c, in, sizeof(cplx) * n * nc);
      for (int j : {0, nc - 1})
        for (int k = 0; k <= n / 2; ++k) {
          const int km = (n - k) % n;
          const cplx h = 0.5 * (in[k * nc + j] + std::conj(in[km * nc + j]));
          c[k * nc + j] = h;
          c[km * nc + j] = std::conj(h);
        }
      fftw_execute_dft_c2r(bwd, cbuf, rbuf);
      std::memcpy(out, rbuf, sizeof(double) * n * n);
    }
  };

  // r2c of both levels; its adjoint; c2r (projected) of both levels; its adjoint
  CVec Rg(const Fft & F, const RVec & x) const {
    const int G = F.n * F.n, M = F.n * F.nc;
    CVec y(2 * M);
    for (int k = 0; k < 2; ++k) F.r2c(x.data() + k * G, y.data() + k * M);
    return y;
  }
  RVec RTg(const Fft & F, const CVec & y) const {  // <r2c x, y> = <x, RT y>
    const int G = F.n * F.n, M = F.n * F.nc;
    CVec z(y);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < F.n; ++i)
        for (int j = 1; j < F.nc - 1; ++j) z[k * M + i * F.nc + j] *= 0.5;
    RVec x(2 * G);
    for (int k = 0; k < 2; ++k) F.c2r(z.data() + k * M, x.data() + k * G);
    return x;
  }
  RVec Cg(const Fft & F, const CVec & y) const {
    const int G = F.n * F.n, M = F.n * F.nc;
    RVec x(2 * G);
    for (int k = 0; k < 2; ++k) F.c2r(y.data() + k * M, x.data() + k * G);
    return x;
  }
  CVec CTg(const Fft & F, const RVec & x) const {  // <c2r y, x> = <y, CT x>
    const int M = F.n * F.nc;
    CVec y = Rg(F, x);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < F.n; ++i)
        for (int j = 1; j < F.nc - 1; ++j) y[k * M + i * F.nc + j] *= 2.;
    return y;
  }
  CVec R(const RVec & x) const { return Rg(fft_, x); }
  RVec RT(const CVec & y) const { return RTg(fft_, y); }
  CVec CT(const RVec & x) const { return CTg(fft_, x); }
  RVec scaleC(const CVec & q) const {  // irfft2 = c2r / N^2
    RVec x = Cg(fft_, q);
    const double norm = 1. / (double(N_) * N_);
    for (double & v : x) v *= norm;
    return x;
  }

  // ------------------------------------------- zero-padding to 3N/2 and truncation
  CVec specpad(const CVec & s) const {
    const int M = N_ * Nc_, Mp = Np_ * Ncp_, h = N_ / 2;
    CVec o(2 * Mp, cplx(0.));
    for (int k = 0; k < 2; ++k) {
      const cplx * src = s.data() + k * M;
      cplx * dst = o.data() + k * Mp;
      for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {
          dst[i * Ncp_ + j] = 2.25 * src[i * Nc_ + j];
          dst[(Np_ - h + i) * Ncp_ + j] = 2.25 * src[(N_ - h + i) * Nc_ + j];
        }
      for (int i = 0; i < h; ++i) {
        dst[i * Ncp_ + h] = std::conj(2.25 * src[i * Nc_ + (Nc_ - 1)]);
        dst[(Np_ - h + i) * Ncp_ + h] = std::conj(2.25 * src[(N_ - h + i) * Nc_ + (Nc_ - 1)]);
      }
    }
    return o;
  }
  CVec specpadT(const CVec & o) const {
    const int M = N_ * Nc_, Mp = Np_ * Ncp_, h = N_ / 2;
    CVec s(2 * M, cplx(0.));
    for (int k = 0; k < 2; ++k) {
      cplx * src = s.data() + k * M;
      const cplx * dst = o.data() + k * Mp;
      for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {
          src[i * Nc_ + j] += 2.25 * dst[i * Ncp_ + j];
          src[(N_ - h + i) * Nc_ + j] += 2.25 * dst[(Np_ - h + i) * Ncp_ + j];
        }
      for (int i = 0; i < h; ++i) {
        src[i * Nc_ + (Nc_ - 1)] += 2.25 * std::conj(dst[i * Ncp_ + h]);
        src[(N_ - h + i) * Nc_ + (Nc_ - 1)] += 2.25 * std::conj(dst[(Np_ - h + i) * Ncp_ + h]);
      }
    }
    return s;
  }
  CVec spectrunc(const CVec & o) const {
    const int M = N_ * Nc_, Mp = Np_ * Ncp_, h = N_ / 2;
    CVec s(2 * M, cplx(0.));
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {
          s[k * M + i * Nc_ + j] = o[k * Mp + i * Ncp_ + j];
          s[k * M + (N_ - h + i) * Nc_ + j] = o[k * Mp + (Np_ - h + i) * Ncp_ + j];
        }
    return s;
  }
  CVec spectruncT(const CVec & s) const {
    const int M = N_ * Nc_, Mp = Np_ * Ncp_, h = N_ / 2;
    CVec o(2 * Mp, cplx(0.));
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < h; ++i)
        for (int j = 0; j < h; ++j) {
          o[k * Mp + i * Ncp_ + j] = s[k * M + i * Nc_ + j];
          o[k * Mp + (Np_ - h + i) * Ncp_ + j] = s[k * M + (N_ - h + i) * Nc_ + j];
        }
    return o;
  }

  // ------------------------------------------------------- PV inversion and adjoint
  CVec invert(const CVec & q) const {
    const int M = N_ * Nc_;
    CVec psi(2 * M);
    for (int i = 0; i < M; ++i) {
      const double h = hovermu_[i], s = sinhmu_[i], t = tanhmu_[i];
      psi[i] = h * (q[M + i] / s - q[i] / t);
      psi[M + i] = h * (q[M + i] / t - q[i] / s);
    }
    return psi;
  }
  CVec invertT(const CVec & pa) const {
    const int M = N_ * Nc_;
    CVec qa(2 * M);
    for (int i = 0; i < M; ++i) {
      const double h = hovermu_[i], s = sinhmu_[i], t = tanhmu_[i];
      qa[i] = -h / t * pa[i] - h / s * pa[M + i];
      qa[M + i] = h / s * pa[i] + h / t * pa[M + i];
    }
    return qa;
  }

  // -------------------------------- x/y derivatives on the padded grid and adjoint
  void deriv(const CVec & q, RVec & dx, RVec & dy) const {
    const int Mp = Np_ * Ncp_;
    const double norm = 1. / (double(Np_) * Np_);
    const CVec p = specpad(q);
    CVec a(2 * Mp), b(2 * Mp);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < Mp; ++i) {
        a[k * Mp + i] = ikp_[i] * p[k * Mp + i];
        b[k * Mp + i] = ilp_[i] * p[k * Mp + i];
      }
    dx = Cg(fftp_, a); dy = Cg(fftp_, b);
    for (double & v : dx) v *= norm;
    for (double & v : dy) v *= norm;
  }
  CVec derivT(const RVec & dxa, const RVec & dya) const {
    const int Mp = Np_ * Ncp_;
    const double norm = 1. / (double(Np_) * Np_);
    RVec xs(dxa), ys(dya);
    for (double & v : xs) v *= norm;
    for (double & v : ys) v *= norm;
    const CVec a = CTg(fftp_, xs), b = CTg(fftp_, ys);
    CVec p(2 * Mp);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < Mp; ++i)
        p[k * Mp + i] = std::conj(ikp_[i]) * a[k * Mp + i] + std::conj(ilp_[i]) * b[k * Mp + i];
    return specpadT(p);
  }

  // ------------------------------------------------------------- tendencies
  // Nonlinear tendency of stage state qs; qn = state at the start of the RK4 step.
  CVec nltend(const CVec & qs, const CVec & qn) const {
    const int M = N_ * Nc_;
    const CVec psi = invert(qs);
    RVec psx, psy, qx, qy;
    deriv(psi, psx, psy); deriv(qs, qx, qy);
    RVec jac(psx.size());
    for (size_t i = 0; i < jac.size(); ++i) jac[i] = psx[i] * qy[i] - psy[i] * qx[i];
    const CVec js = spectrunc(Rg(fftp_, jac));
    CVec f(2 * M);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < M; ++i) {
        const int I = k * M + i;
        f[I] = (pveq_[I] - qs[I]) / tdiab_ - js[I] + r_[k] * ksq_[i] * psi[I]
               + hyperdiff_[i] * qn[I];
      }
    return f;
  }
  // TL tendency: dqs = stage perturbation, dqn = step-start perturbation, qs = stage state
  CVec tltend(const CVec & dqs, const CVec & dqn, const CVec & qs) const {
    const int M = N_ * Nc_;
    RVec psx, psy, qx, qy, dpsx, dpsy, dqx, dqy;
    deriv(invert(qs), psx, psy); deriv(qs, qx, qy);
    const CVec dpsi = invert(dqs);
    deriv(dpsi, dpsx, dpsy); deriv(dqs, dqx, dqy);
    RVec djac(psx.size());
    for (size_t i = 0; i < djac.size(); ++i)
      djac[i] = dpsx[i] * qy[i] + psx[i] * dqy[i] - dpsy[i] * qx[i] - psy[i] * dqx[i];
    const CVec djs = spectrunc(Rg(fftp_, djac));
    CVec df(2 * M);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < M; ++i) {
        const int I = k * M + i;
        df[I] = -dqs[I] / tdiab_ - djs[I] + r_[k] * ksq_[i] * dpsi[I] + hyperdiff_[i] * dqn[I];
      }
    return df;
  }
  // Adjoint of tltend: dfa -> (dqsa, dqna)
  void adtend(const CVec & dfa, const CVec & qs, CVec & dqsa, CVec & dqna) const {
    const int M = N_ * Nc_;
    RVec psx, psy, qx, qy;
    deriv(invert(qs), psx, psy); deriv(qs, qx, qy);
    dqsa.assign(2 * M, cplx(0.)); dqna.assign(2 * M, cplx(0.));
    CVec dpsia(2 * M), djsa(2 * M);
    for (int k = 0; k < 2; ++k)
      for (int i = 0; i < M; ++i) {
        const int I = k * M + i;
        dqsa[I] = -dfa[I] / tdiab_;
        dpsia[I] = r_[k] * ksq_[i] * dfa[I];
        dqna[I] = hyperdiff_[i] * dfa[I];
        djsa[I] = -dfa[I];
      }
    const RVec djaca = RTg(fftp_, spectruncT(djsa));
    RVec dpsxa(djaca.size()), dpsya(djaca.size()), dqxa(djaca.size()), dqya(djaca.size());
    for (size_t i = 0; i < djaca.size(); ++i) {
      dpsxa[i] = djaca[i] * qy[i];
      dqya[i] = djaca[i] * psx[i];
      dpsya[i] = -djaca[i] * qx[i];
      dqxa[i] = -djaca[i] * psy[i];
    }
    add(dpsia, derivT(dpsxa, dpsya));
    add(dqsa, derivT(dqxa, dqya));
    add(dqsa, invertT(dpsia));
  }

  // One RK4 step of the nonlinear model; optionally keeps the four stage states.
  CVec rk4(const CVec & q, CVec * stages) const {
    const double a = 0.5 * dt_;
    const CVec k1 = nltend(q, q);
    const CVec q2 = axpy(q, a, k1);
    const CVec k2 = nltend(q2, q);
    const CVec q3 = axpy(q, a, k2);
    const CVec k3 = nltend(q3, q);
    const CVec q4 = axpy(q, dt_, k3);
    const CVec k4 = nltend(q4, q);
    if (stages) { stages[0] = q; stages[1] = q2; stages[2] = q3; stages[3] = q4; }
    CVec out(q);
    for (size_t i = 0; i < out.size(); ++i)
      out[i] += dt_ / 6. * (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]);
    return out;
  }

  static CVec axpy(const CVec & x, double a, const CVec & y) {
    CVec z(x);
    for (size_t i = 0; i < z.size(); ++i) z[i] += a * y[i];
    return z;
  }
  static CVec scale(const CVec & x, double a) {
    CVec z(x);
    for (cplx & v : z) v *= a;
    return z;
  }
  static void add(CVec & x, const CVec & y) {
    for (size_t i = 0; i < x.size(); ++i) x[i] += y[i];
  }

  int N_, Nc_, Np_, Ncp_;
  double tdiab_, dt_, r_[2];
  std::vector<double> ksq_, hovermu_, tanhmu_, sinhmu_, hyperdiff_;
  CVec ikp_, ilp_, pveq_;
  Fft fft_, fftp_;
};
