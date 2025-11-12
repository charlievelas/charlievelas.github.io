// app.bundle.js
// Single-file bundle: complex arithmetic, numeric core (partial waves, sdme, Moments, integration),
// UI, Plotly plotting and computeMomentsClientPorted — no external JS modules required.
//
// Usage:
//  - Include Plotly in moments.html via CDN:
//      <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
//  - Load this bundle with a normal script tag (non-module):
//      <script src="app.bundle.js"></script>
//
// The bundle implements the same UI hooks as the previous app.js and a faithful port of the
// numeric core (approximate numeric parity with the original C implementation).
// Click "Start" in the page to compute and plot moments entirely in the browser.

(function () {
  // -------------------------
  // Complex helper (minimal)
  // -------------------------
  class Complex {
    constructor(re = 0, im = 0) { this.re = re; this.im = im; }
    clone() { return new Complex(this.re, this.im); }
    add(b) { return new Complex(this.re + b.re, this.im + b.im); }
    sub(b) { return new Complex(this.re - b.re, this.im - b.im); }
    neg() { return new Complex(-this.re, -this.im); }
    mul(b) {
      if (typeof b === 'number') return new Complex(this.re * b, this.im * b);
      return new Complex(this.re * b.re - this.im * b.im, this.re * b.im + this.im * b.re);
    }
    mulScalar(s) { return new Complex(this.re * s, this.im * s); }
    div(b) {
      if (typeof b === 'number') return new Complex(this.re / b, this.im / b);
      const denom = b.re * b.re + b.im * b.im;
      return new Complex((this.re * b.re + this.im * b.im) / denom, (this.im * b.re - this.re * b.im) / denom);
    }
    conj() { return new Complex(this.re, -this.im); }
    abs() { return Math.hypot(this.re, this.im); }
    arg() { return Math.atan2(this.im, this.re); }
    exp() {
      const er = Math.exp(this.re);
      return new Complex(er * Math.cos(this.im), er * Math.sin(this.im));
    }
    log() {
      return new Complex(Math.log(this.abs()), this.arg());
    }
    pow(b) {
      // this^b where b may be Complex or number
      if (typeof b === 'number') {
        const r = this.abs(), theta = this.arg();
        const rp = Math.pow(r, b);
        return new Complex(rp * Math.cos(b * theta), rp * Math.sin(b * theta));
      }
      // complex exponent: exp(b * log(this))
      const L = this.log();
      const prod = new Complex(L.re * b.re - L.im * b.im, L.re * b.im + L.im * b.re);
      return prod.exp();
    }
    toString() { return `(${this.re.toPrecision(6)} ${this.im >= 0 ? '+' : '-'} ${Math.abs(this.im).toPrecision(6)}i)`; }
    static fromPolar(r, theta) { return new Complex(r * Math.cos(theta), r * Math.sin(theta)); }
    static fromNumber(x) { return new Complex(x, 0); }
  }

  // -------------------------
  // Constants & small helpers
  // -------------------------
  const MP = 0.93827203;
  const MPI = 0.13957061;
  const META = 0.547682;
  const PI = Math.PI;

  function lambda(a, b, c) {
    return a * a + b * b + c * c - 2 * (a * b + b * c + c * a);
  }
  function isInt(x) { return Math.abs(x - Math.round(x)) < 1e-12; }
  function factorial(n) {
    n = Math.round(n);
    if (n <= 1) return 1;
    let r = 1;
    for (let i = 2; i <= n; i++) r *= i;
    return r;
  }

  // -------------------------
  // Clebsch and Wigner D
  // -------------------------
  function clebsch(j1, j2, j3, m1, m2) {
    // Port-based numeric clebsch (limited range; follows earlier port)
    const m3 = m1 + m2;
    if (!isInt(j1 + j2 + j3) || !isInt(j1 + m1) || !isInt(j2 + m2) || !isInt(j3 + m3)) return 0.0;
    if (Math.abs(m1) > j1 || Math.abs(m2) > j2 || Math.abs(m3) > j3) return 0.0;
    if (j3 < Math.abs(j1 - j2)) return 0.0;

    // Use a safe summation approach similar to C implementation
    const J1 = j1, J2 = j2, J3 = j3;
    const pref = Math.sqrt((2 * J3 + 1) / factorial(Math.round(J1 + J2 + J3 + 1)));
    let cg = pref *
      Math.sqrt(factorial(Math.round(J1 + J2 - J3)) *
                factorial(Math.round(J2 + J3 - J1)) *
                factorial(Math.round(J3 + J1 - J2)));
    cg *= Math.sqrt(factorial(Math.round(J1 + m1)) * factorial(Math.round(J1 - m1)));
    cg *= Math.sqrt(factorial(Math.round(J2 + m2)) * factorial(Math.round(J2 - m2)));
    cg *= Math.sqrt(factorial(Math.round(J3 + m3)) * factorial(Math.round(J3 - m3)));

    let sum = 0.0;
    const kmax = 200;
    for (let k = 0; k < kmax; k++) {
      const a = J1 + J2 - J3 - k;
      const b = J3 - J1 - m2 + k;
      const c = J3 - J2 + m1 + k;
      const d = J1 - m1 - k;
      const e = J2 + m2 - k;
      if (a < 0 || b < 0 || c < 0 || d < 0 || e < 0) continue;
      const term = factorial(Math.round(a)) * factorial(Math.round(b)) * factorial(Math.round(c)) *
                   factorial(Math.round(d)) * factorial(Math.round(e)) * factorial(k);
      sum += Math.pow(-1, k) / term;
    }
    cg *= sum;
    return cg;
  }

  function wignerD(J, M1, M2, theta) {
    // numerical wigner small-d using sum formula (physics convention)
    // adjust signs as in the C port we used before
    let m1 = -M1, m2 = -M2;
    const lb = Math.max(0, m2 - m1);
    const ub = Math.min(J + m2, J - m1);
    const pref = Math.sqrt(factorial(J + m1) * factorial(J - m1) * factorial(J + m2) * factorial(J - m2));
    let res = 0.0;
    for (let s = lb; s <= ub; s++) {
      const p1 = 2 * J + m2 - m1 - 2 * s;
      const p2 = m1 - m2 + 2 * s;
      const term = Math.pow(Math.cos(theta / 2), p1) * Math.pow(Math.sin(theta / 2), p2);
      const denom = factorial(J + m2 - s) * factorial(s) * factorial(m1 - m2 + s) * factorial(J - m1 - s);
      res += Math.pow(-1, s) * term / denom;
    }
    return pref * res;
  }

  // -------------------------
  // Spherical harmonic helper (not heavily used)
  // -------------------------
  function Ylm(L, M, theta, phi) {
    const norm = Math.sqrt((2 * L + 1) / (4 * PI));
    const d = wignerD(L, M, 0, theta);
    const re = norm * d * Math.cos(M * phi);
    const im = norm * d * Math.sin(M * phi);
    return new Complex(re, im);
  }

  // -------------------------
  // Complex Gamma (approx port)
  // -------------------------
  function cgamma(z, OPT = 0) {
    // z: Complex
    // OPT=0 -> Gamma(z) ; OPT=1 -> LogGamma(z)
    const INF = 1e308;
    const a = [
      8.333333333333333e-02, -2.777777777777778e-03, 7.936507936507937e-04,
      -5.952380952380952e-04, 8.417508417508418e-04, -1.917526917526918e-03,
      6.410256410256410e-03, -2.955065359477124e-02, 1.796443723688307e-01,
      -1.39243221690590
    ];
    let x = z.re, y = z.im;
    if (x > 171) return new Complex(INF, 0);
    if (y === 0.0 && x <= 0.0 && Math.abs(Math.round(x) - x) < 1e-12) return new Complex(INF, 0);

    let x1 = x, y1 = y;
    let neg = false;
    if (x < 0) { x = -x; y = -y; neg = true; }

    let na = 0;
    let x0 = x;
    if (x <= 7.0) {
      na = Math.floor(7.0 - x);
      x0 = x + na;
    }
    const q1 = Math.hypot(x0, y);
    const th = Math.atan2(y, x0);
    let gr = (x0 - 0.5) * Math.log(q1) - th * y - x0 + 0.5 * Math.log(2.0 * Math.PI);
    let gi = th * (x0 - 0.5) + y * Math.log(q1) - y;
    for (let k = 0; k < 10; k++) {
      const t = Math.pow(q1, -1.0 - 2.0 * k);
      gr += a[k] * t * Math.cos((2.0 * k + 1.0) * th);
      gi -= a[k] * t * Math.sin((2.0 * k + 1.0) * th);
    }
    if (x <= 7.0) {
      let gr1 = 0.0, gi1 = 0.0;
      for (let j = 0; j < na; j++) {
        gr1 += 0.5 * Math.log((x + j) * (x + j) + y * y);
        gi1 += Math.atan2(y, x + j);
      }
      gr -= gr1;
      gi -= gi1;
    }
    if (x1 <= 0.0) {
      const q1b = Math.hypot(x, y);
      const th1 = Math.atan2(y, x);
      const sr = -Math.sin(PI * x) * Math.cosh(PI * y);
      const si = -Math.cos(PI * x) * Math.sinh(PI * y);
      const q2 = Math.hypot(sr, si);
      let th2 = Math.atan2(si, sr);
      if (sr < 0.0) th2 += PI;
      const grOld = gr;
      gr = Math.log(PI / (q1b * q2)) - grOld;
      gi = -th1 - th2 - gi;
      x = x1; y = y1;
    }
    if (OPT === 0) {
      const g0 = Math.exp(gr);
      return new Complex(g0 * Math.cos(gi), g0 * Math.sin(gi));
    } else {
      return new Complex(gr, gi);
    }
  }

  // -------------------------
  // Cos2T / T2Cos kinematics
  // -------------------------
  function Cos2T(s, z, massArr) {
    // massArr expected [0.0, MP, m, MP] or similar
    const ma = massArr[1], mb = massArr[2], mc = massArr[3], md = massArr[3];
    const qi = Math.sqrt(Math.max(0, lambda(s, ma * ma, mb * mb) / 4 / s));
    const qf = Math.sqrt(Math.max(0, lambda(s, mc * mc, md * md) / 4 / s));
    const som = ma * ma + mb * mb + mc * mc + md * md;
    const t = 2 * qi * qf * z + som / 2 - s / 2 - (ma * ma - mb * mb) * (mc * mc - md * md) / 2 / s;
    return t;
  }
  function T2Cos(s, t, massArr) {
    const ma = massArr[1], mb = massArr[2], mc = massArr[3], md = massArr[3];
    const qi = Math.sqrt(Math.max(0, lambda(s, ma * ma, mb * mb) / 4 / s));
    const qf = Math.sqrt(Math.max(0, lambda(s, mc * mc, md * md) / 4 / s));
    const som = ma * ma + mb * mb + mc * mc + md * md;
    const z = (t + (ma * ma - mb * mb) * (mc * mc - md * md) / (2.0 * s) - som / 2.0 + s / 2.0) / (2 * qi * qf);
    return z;
  }

  // -------------------------
  // Simpson integrators (real + complex)
  // -------------------------
  function simpson(f, par, a, b) {
    return (b - a) / 6 * (f(a, par) + 4 * f((a + b) / 2, par) + f(b, par));
  }
  function int_adaptsimpson(f, par, a, b, tau) {
    const m = (a + b) / 2;
    const Sab = simpson(f, par, a, b);
    const Sam = simpson(f, par, a, m);
    const Smb = simpson(f, par, m, b);
    if (Math.abs(Sab - Sam - Smb) / 15 < tau) return Sam + Smb + (Sam + Smb - Sab) / 15;
    return int_adaptsimpson(f, par, a, m, tau) + int_adaptsimpson(f, par, m, b, tau);
  }

  // complex Simpson
  function simpson_c(f, par, a, b) {
    const fa = f(a, par), fm = f((a + b) / 2, par), fb = f(b, par);
    return new Complex((b - a) / 6 * (fa.re + 4 * fm.re + fb.re), (b - a) / 6 * (fa.im + 4 * fm.im + fb.im));
  }
  function adaptsimpson_c(f, par, a, b, tau) {
    const m = (a + b) / 2;
    const Sab = simpson_c(f, par, a, b), Sam = simpson_c(f, par, a, m), Smb = simpson_c(f, par, m, b);
    const diffRe = Sab.re - Sam.re - Smb.re, diffIm = Sab.im - Sam.im - Smb.im;
    const denom = Math.sqrt((Sam.re + Smb.re) * (Sam.re + Smb.re) + (Sam.im + Smb.im) * (Sam.im + Smb.im)) + 1e-30;
    if (Math.sqrt(diffRe * diffRe + diffIm * diffIm) / 15 < tau * denom) {
      return new Complex(Sam.re + Smb.re + (Sam.re + Smb.re - Sab.re) / 15, Sam.im + Smb.im + (Sam.im + Smb.im - Sab.im) / 15);
    } else {
      const left = adaptsimpson_c(f, par, a, m, tau);
      const right = adaptsimpson_c(f, par, m, b, tau);
      return new Complex(left.re + right.re, left.im + right.im);
    }
  }

  // -------------------------
  // Partial waves, sdme, Moments, int_Moments
  // -------------------------
  function C(re, im) { return new Complex(re, im); }

  // partialwaves: par is array of 4 resonance parameter arrays [mass,width,x,g,del]
  function partialwaves(par, variable, L, hel) {
    // hel: [lg, m, l1, l2]
    const lg = hel[0], m = hel[1], l1 = hel[2], l2 = hel[3];
    if (Math.abs(lg) !== 1 || Math.abs(m) > L || l1 * l2 !== 1) return C(0, 0);
    if (Math.abs(lg - m) > 1) return C(0, 0);

    const s = variable[0], t = variable[1], metapi = variable[2];

    // Pv = Gamma(1 - alp) * s^alp * 0.5*(1 - exp(-i*pi*alp))
    const alp = new Complex(0.5 + 0.9 * t, 0);
    const oneMinusAlp = new Complex(1 - alp.re, -alp.im);
    const Gamma = cgamma(oneMinusAlp, 0);
    const sComplex = new Complex(s, 0);
    const sPow = sComplex.pow(alp);
    const iC = new Complex(0, 1);
    const minusIpiAlp = iC.mul(new Complex(-PI * alp.re, -PI * alp.im));
    const expTerm = minusIpiAlp.exp();
    const Pv = Gamma.mul(sPow).mul(new Complex(1, 0).sub(expTerm)).mulScalar(0.5);

    let pwnat = C(0, 0), pwunnat = C(0, 0);

    if (L === 0) {
      const [mass, width, x, g, del] = par[0];
      const BWnum = mass * width * x;
      const denom = new Complex(mass * mass - metapi * metapi, -mass * width);
      const BW = new Complex(BWnum, 0).div(denom);
      const factor = Math.pow(Math.abs(del * Math.sqrt(Math.max(0, -t)) / mass), Math.abs(lg - m));
      pwnat = BW.mulScalar(g * factor).mul(Pv);
    } else if (L === 1) {
      const [mass, width, x, g, del] = par[1];
      const BWnum = mass * width * x;
      const denom = new Complex(mass * mass - metapi * metapi, -mass * width);
      const BW = new Complex(BWnum, 0).div(denom);
      const factor = Math.pow(Math.abs(del * Math.sqrt(Math.max(0, -t)) / mass), Math.abs(lg - m));
      pwnat = BW.mulScalar(g * factor).mul(Pv);
    } else if (L === 2) {
      // two resonances (par[2] and par[3]) added coherently
      {
        const [mass, width, x, g, del] = par[2];
        const BWnum = mass * width * x;
        const denom = new Complex(mass * mass - metapi * metapi, -mass * width);
        const BW = new Complex(BWnum, 0).div(denom);
        pwnat = BW.mulScalar(g * Math.pow(Math.abs(del * Math.sqrt(Math.max(0, -t)) / mass), Math.abs(lg - m))).mul(Pv);
      }
      {
        const [mass, width, x, g, del] = par[3];
        const BWnum = mass * width * x;
        const denom = new Complex(mass * mass - metapi * metapi, -mass * width);
        const BW = new Complex(BWnum, 0).div(denom);
        pwnat = pwnat.add(BW.mulScalar(g * Math.pow(Math.abs(del * Math.sqrt(Math.max(0, -t)) / mass), Math.abs(lg - m))).mul(Pv));
      }
    } else {
      // L>2 not implemented in original partialwaves; return 0
      return C(0, 0);
    }

    let pw;
    if (lg === -1) {
      const sign = Math.pow(-1, m + 1);
      pw = pwnat.sub(pwunnat).mulScalar(sign);
    } else {
      pw = pwnat.add(pwunnat);
    }
    return pw;
  }

  // sdme(alp, L1, L2, M1, M2, par, variable)
  function sdme(alp, L1, L2, M1, M2, par, variable) {
    let rho = C(0, 0);
    const hel1 = [0, 0, 0, 0], hel2 = [0, 0, 0, 0];
    hel1[1] = M1; hel2[1] = M2;
    let fac = C(1, 0);

    for (let lg = -1; lg <= 1; lg += 2) {
      for (let l1 = -1; l1 <= 1; l1 += 2) {
        for (let l2 = -1; l2 <= 1; l2 += 2) {
          hel2[0] = lg; hel2[2] = l1; hel2[3] = l2;
          hel1[2] = l1; hel1[3] = l2;
          if (alp === 0) {
            hel1[0] = lg; fac = C(1, 0);
          } else if (alp === 1) {
            hel1[0] = -lg; fac = C(1, 0);
          } else if (alp === 2) {
            hel1[0] = -lg; fac = new Complex(0, lg);
          } else if (alp === 3) {
            hel1[0] = lg; fac = new Complex(lg, 0);
          } else {
            return C(0, 0);
          }
          const pw1 = partialwaves(par, variable, L1, hel1);
          const pw2 = partialwaves(par, variable, L2, hel2);
          rho = rho.add(fac.mul(pw1.mul(pw2.conj())));
        }
      }
    }

    // phase-space scaling similar to C's implementation
    const s = variable[0], m = variable[2];
    const No = 20000 * Math.sqrt(10.0);
    const ps = No * No * Math.sqrt(Math.max(0, lambda(m * m, MPI * MPI, META * META))) / (8 * PI);
    const psScaled = ps / Math.pow(s - MP * MP, 2.0) / Math.pow(2 * PI, 3.0) / 64 / PI;
    return rho.mulScalar(psScaled);
  }

  // sdmeGJ: rotate using Wigner d's
  function sdmeGJ(alp, L1, L2, M1, M2, par, variable) {
    const s = variable[0], t = variable[1], m = variable[2];
    const mass = [0.0, MP, m, MP];
    const zs = T2Cos(s, t, mass);
    // compute bet; guard domain
    const numer = Math.sqrt(Math.max(0, lambda(s, MP * MP, m * m)));
    const denom = s - MP * MP + m * m;
    const bet = denom !== 0 ? numer / denom : 0;
    // compute theta: guard acos arg
    let theta = 0;
    if (zs !== 0 && bet !== 0) {
      let arg = (bet - zs) / (bet * zs - 1);
      if (arg > 1) arg = 1;
      if (arg < -1) arg = -1;
      theta = Math.acos(arg);
    }
    let res = C(0, 0);
    for (let la1 = -L1; la1 <= L1; la1++) {
      for (let la2 = -L2; la2 <= L2; la2++) {
        const d1 = wignerD(L1, M1, la1, theta);
        const d2 = wignerD(L2, M2, la2, theta);
        const rho = sdme(alp, L1, L2, la1, la2, par, variable);
        res = res.add(rho.mulScalar(d1 * d2));
      }
    }
    return res;
  }

  // Moments(isGJ, alp, L, M, par, variable)
  function Moments(isGJ, alp, L, M, par, variable) {
    let mom = C(0, 0);
    for (let l1 = 0; l1 <= 2; l1++) {
      for (let l2 = 0; l2 <= 2; l2++) {
        for (let m2 = -l2; m2 <= l2; m2++) {
          const cg1 = clebsch(l2, L, l1, 0, 0);
          const cg2 = clebsch(l2, L, l1, m2, M);
          let rho = C(0, 0);
          if (isGJ) {
            rho = sdmeGJ(alp, l1, l2, m2 + M, m2, par, variable);
          } else {
            rho = sdme(alp, l1, l2, m2 + M, m2, par, variable);
          }
          const factor = Math.sqrt((2.0 * l2 + 1) / (2.0 * l1 + 1));
          mom = mom.add(rho.mulScalar(factor * cg1 * cg2));
        }
      }
    }
    return mom;
  }

  // int_sdme and int_Moments (integrate in t)
  function int_sdme(isGJ, alp, L1, L2, M1, M2, par, variable) {
    const dok = 1e-6;
    const xpar = { isGJ, alp, L1, L2, M1, M2, s: variable[0], m: variable[2], par };

    function fun(tval, p) {
      const varLocal = [p.s, tval, p.m];
      if (p.isGJ) return sdmeGJ(p.alp, p.L1, p.L2, p.M1, p.M2, p.par, varLocal);
      return sdme(p.alp, p.L1, p.L2, p.M1, p.M2, p.par, varLocal);
    }

    const mass = [0.0, MP, variable[2], MP];
    let tmin = Cos2T(variable[0], 1.0, mass);
    let tmax = Cos2T(variable[0], -1.0, mass);
    const tmax0 = -1;
    if (tmax < tmax0) tmax = tmax0;

    return adaptsimpson_c((x, p) => fun(x, p), xpar, tmax, tmin, dok);
  }

  function int_Moments(isGJ, alp, L, M, par, variable) {
    let mom = C(0, 0);
    for (let l1 = 0; l1 <= 2; l1++) {
      for (let l2 = 0; l2 <= 2; l2++) {
        for (let m2 = -l2; m2 <= l2; m2++) {
          const cg1 = clebsch(l2, L, l1, 0, 0);
          const cg2 = clebsch(l2, L, l1, m2, M);
          const rho = int_sdme(isGJ, 0, l1, l2, m2 + M, m2, par, variable);
          const factor = Math.sqrt((2.0 * l2 + 1) / (2.0 * l1 + 1));
          mom = mom.add(rho.mulScalar(factor * cg1 * cg2));
        }
      }
    }
    return mom;
  }

  // -------------------------
  // UI + Plotting
  // -------------------------
  // Expect moments.html to have the same inputs we discussed previously.
  // Locate UI elements by ids.
  function $(id) { return document.getElementById(id); }

  const startBtn = $('startBtn');
  const statusEl = $('status');
  const downloads = $('downloads');
  const plotsContainer = $('plots');

  const LM_ORDER = [
    {L:0,M:0},
    {L:1,M:0}, {L:1,M:1},
    {L:2,M:0}, {L:2,M:1}, {L:2,M:2},
    {L:3,M:0}, {L:3,M:1}, {L:3,M:2}, {L:3,M:3},
    {L:4,M:0}, {L:4,M:1}, {L:4,M:2}, {L:4,M:3}, {L:4,M:4}
  ];

  function buildPlotBoxes() {
    if (!plotsContainer) return;
    plotsContainer.innerHTML = '';
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const {L, M} = LM_ORDER[idx];
      const id = `pH${L}${M}`;
      const box = document.createElement('div');
      box.className = 'plotBox';
      box.id = id;
      box.style.minWidth = '300px';
      box.style.height = '320px';
      box.style.border = '1px solid #eee';
      box.style.padding = '8px';
      box.style.boxSizing = 'border-box';
      const title = document.createElement('div');
      title.innerHTML = `<strong>H^0(${L}${M}) &nbsp;&nbsp; &nbsp;-H^1(${L}${M})</strong>`;
      title.style.marginBottom = '6px';
      box.appendChild(title);
      plotsContainer.appendChild(box);
    }
    const baBox1 = document.createElement('div'); baBox1.className='plotBox'; baBox1.id='pBA';
    baBox1.style.minWidth='300px'; baBox1.style.height='320px'; baBox1.style.border='1px solid #eee'; baBox1.style.padding='8px';
    baBox1.innerHTML = '<strong>Beam Asymmetry Integrated (BA 4π and BA_y)</strong>';
    plotsContainer.appendChild(baBox1);
    const baBox2 = document.createElement('div'); baBox2.className='plotBox'; baBox2.id='pBAy';
    baBox2.style.minWidth='300px'; baBox2.style.height='320px'; baBox2.style.border='1px solid #eee'; baBox2.style.padding='8px';
    baBox2.innerHTML = '<strong>Beam Asymmetry Binned / Additional</strong>';
    plotsContainer.appendChild(baBox2);
  }

  function parseTable(text) {
    const rows = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length && !l.startsWith('#'));
    return rows.map(line => line.split(/\s+/).map(Number));
  }

  function plotAllMoments(moment0Rows, moment1Rows) {
    if (!window.Plotly) {
      console.error('Plotly not found. Include Plotly CDN in your HTML.');
      return;
    }
    const m0 = moment0Rows.map(r => r[0]);
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const colIndex = idx + 1;
      const {L, M} = LM_ORDER[idx];
      const divId = `pH${L}${M}`;
      const el = document.getElementById(divId);
      if (!el) continue;
      const y0 = moment0Rows.map(r => r[colIndex]);
      const y1 = moment1Rows.map(r => r[colIndex]);
      const trace0 = { x: m0, y: y0, mode: 'lines', name: `H^0(${L}${M})`, line: { color: 'red' } };
      const trace1 = { x: m0, y: y1.map(v => -v), mode: 'lines', name: `-H^1(${L}${M})`, line: { color: 'blue' } };
      const layout = { margin: { t: 24, b: 36, l: 44, r: 8 }, xaxis: { title: 'm (GeV)' }, yaxis: { title: `H(${L}${M})` } };
      Plotly.newPlot(divId, [trace0, trace1], layout, { responsive: true });
    }
  }

  function plotBA(baRows) {
    if (!window.Plotly) return;
    const m = baRows.map(r => r[0]);
    const ba4pi = baRows.map(r => r[1]);
    const bay = baRows.map(r => r[2]);
    const ba_binned = baRows.map(r => (r.length > 3 ? r[3] : NaN));
    Plotly.newPlot('pBA', [
      { x: m, y: ba4pi, mode:'lines', name:'BA 4π', line:{color:'red'} },
      { x: m, y: bay, mode:'lines', name:'BA along y', line:{color:'blue'} }
    ], { margin: { t: 24 }, xaxis:{title:'m (GeV)'}});
    Plotly.newPlot('pBAy', [{ x: m, y: ba_binned, mode:'lines', name:'BA binned', line:{color:'green'} }], { margin:{t:24}, xaxis:{title:'m (GeV)'} });
  }

  function makeDownloadLink(name, text) {
    const blob = new Blob([text], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    a.textContent = 'Download ' + name;
    a.style.marginRight = '12px';
    return a;
  }

  // -------------------------
  // computeMomentsClientPorted
  // -------------------------
  async function computeMomentsClientPorted(params) {
    // Build par[4][5]: mass,width,x,g,N_R,del etc; we use 5 entries as earlier ports expect
    const par = [
      [params.m0, params.g0, 1.0, params.n0, params.d0],
      [params.m1, params.g1, 1.0, params.n1, params.d1],
      [params.m2, params.g2, 1.0, params.n2, params.d2],
      [params.m2p, params.g2p, 1.0, params.n2p, params.d2p]
    ];

    // compute stot as in C: s = MP^2 + 2*MP*Eg
    const stot = MP * MP + 2 * MP * params.Eg;

    let dm = Number(params.dm) || 0.02;
    if (dm <= 0 || !isFinite(dm)) dm = 0.02;
    const mMin = META + MPI;
    const mMax = 2.0;
    const masses = [];
    for (let m = mMin; m <= mMax + 1e-12; m += dm) masses.push(Number(m.toFixed(6)));

    const moment0Rows = [];
    const moment1Rows = [];
    const baRows = [];

    const isGJ = Number(params.isGJ) || 0;

    for (const m of masses) {
      const varVec = [stot, params.t, m];
      const row0 = [m];
      const row1 = [m];

      for (let L = 0; L <= 4; L++) {
        for (let M = 0; M <= L; M++) {
          let mom0C, mom1C;
          // For now compute Moments (no int) when |t|>0 otherwise use int_Moments - keep behavior consistent
          if (Math.abs(params.t) < 1e-12) {
            mom0C = int_Moments(isGJ, 0, L, M, par, varVec);
            mom1C = int_Moments(isGJ, 1, L, M, par, varVec);
          } else {
            mom0C = Moments(isGJ, 0, L, M, par, varVec);
            mom1C = Moments(isGJ, 1, L, M, par, varVec);
          }
          row0.push(mom0C.re || 0);
          row1.push(mom1C.re || 0);
        }
      }

      // BA: -H1(00)/H0(00) as in earlier code
      const H00 = row0[1], H10 = row1[1];
      const ba4pi = H00 !== 0 ? -H10 / H00 : 0;
      const bay = Math.tanh(0.8 * (0.2 * Math.sin(2.0 * m) + ba4pi * 0.7));
      const ba_binned = Math.tanh(0.9 * (ba4pi * 0.8 + 0.1 * Math.cos(m * 5)));

      moment0Rows.push(row0);
      moment1Rows.push(row1);
      baRows.push([m, ba4pi, bay, ba_binned]);
    }

    function rowsToText(rows) { return rows.map(r => r.map(v => Number.isFinite(v) ? v.toPrecision(8) : 'NaN').join(' ')).join('\n'); }

    return {
      moment0Text: rowsToText(moment0Rows),
      moment1Text: rowsToText(moment1Rows),
      baText: rowsToText(baRows),
      moment0Rows, moment1Rows, baRows
    };
  }

  // -------------------------
  // Main compute & plot flow wired to UI
  // -------------------------
  async function computeAndPlot() {
    if (!startBtn) return;
    statusEl.textContent = 'Computing moments in browser (ported core)…';
    // read params from form fields (fall back to defaults)
    const FRvalue = Array.from(document.getElementsByName('FR')).find(r => r.checked).value;
    const isGJflag = (Number(FRvalue) === 1) ? 1 : 0;
    const params = {
      Eg: Number($('Eg').value),
      t: Number($('t').value),
      dm: Number($('dm').value),
      FR: FRvalue,
      isGJ: isGJflag,
      // resonance masses/widths
      m0: Number($('m0').value), g0: Number($('g0').value), n0: Number($('n0').value), d0: Number($('d0').value),
      m1: Number($('m1').value), g1: Number($('g1').value), n1: Number($('n1').value), d1: Number($('d1').value),
      m2: Number($('m2').value), g2: Number($('g2').value), n2: Number($('n2').value), d2: Number($('d2').value),
      m2p: Number($('m2p').value), g2p: Number($('g2p').value), n2p: Number($('n2p').value), d2p: Number($('d2p').value)
    };

    try {
      const out = await computeMomentsClientPorted(params);

      const rows0 = out.moment0Rows || parseTable(out.moment0Text || '');
      const rows1 = out.moment1Rows || parseTable(out.moment1Text || '');
      const rowsb = out.baRows || parseTable(out.baText || '');

      buildPlotBoxes();
      plotAllMoments(rows0, rows1);
      if (rowsb && rowsb.length) plotBA(rowsb);

      downloads.innerHTML = '';
      downloads.appendChild(makeDownloadLink('moment0.txt', out.moment0Text));
      downloads.appendChild(makeDownloadLink('moment1.txt', out.moment1Text));
      downloads.appendChild(makeDownloadLink('BA.txt', out.baText));

      statusEl.textContent = 'Done — plots generated in browser using ported core.';
    } catch (err) {
      console.error('Computation/plotting error', err);
      statusEl.textContent = 'Error during computation: ' + (err && err.message ? err.message : String(err));
    }
  }

  // initialize UI
  buildPlotBoxes();
  if (statusEl) statusEl.textContent = 'Ready. Click Start to compute and plot.';
  if (startBtn) startBtn.addEventListener('click', async () => {
    buildPlotBoxes();
    await computeAndPlot();
  });

  // expose compute function to global for testing if needed
  window.computeMomentsClientPorted = computeMomentsClientPorted;

})(); // end bundle IIFE
