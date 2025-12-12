// bw_ported.js
// Port of the C++ Breit-Wigner + wave_builder + all_moments_builder
// Exposes: window.bwFunc, window.waveBuilder, window.allMomentsBuilder

(function(){
  // Small complex number helper
  class C {
    constructor(re=0, im=0){ this.re = Number(re) || 0; this.im = Number(im) || 0; }
    clone(){ return new C(this.re, this.im); }
    add(b){ return new C(this.re + b.re, this.im + b.im); }
    mul(b){
      if (typeof b === 'number') return new C(this.re * b, this.im * b);
      return new C(this.re*b.re - this.im*b.im, this.re*b.im + this.im*b.re);
    }
    div(b){
      if (typeof b === 'number') return new C(this.re / b, this.im / b);
      const denom = b.re*b.re + b.im*b.im || 1e-30;
      return new C((this.re*b.re + this.im*b.im)/denom, (this.im*b.re - this.re*b.im)/denom);
    }
    conj(){ return new C(this.re, -this.im); }
  }

  function norm(z){ return (z.re*z.re + z.im*z.im); }
  function realProd(z1, z2){ // real( z1 * conj(z2) )
    return z1.re * z2.re + z1.im * z2.im;
  }

  // Breit-Wigner function (ported)
  // massRes, widthRes, massKK, x (amplitude multiplier)
  function bwFunc(massRes, widthRes, massKK, x){
    massRes = Number(massRes);
    widthRes = Number(widthRes);
    massKK = Number(massKK);
    x = (typeof x === 'undefined') ? 1.0 : Number(x);
    const num = new C(massRes * widthRes * x, 0);
    const denom = new C(massRes*massRes - massKK*massKK, - massRes * widthRes);
    return num.div(denom);
  }

  // waveBuilder: resonances is an array of objects with fields {L, m, mass, width, amplitude}
  // massKK is the running mass (from HTML input); if a resonance object has an amplitude field it will be used
  function waveBuilder(resonances, massKK, xOverride){
    massKK = Number(massKK);
    const S0 = new C(0,0);
    let P0 = new C(0,0), P1 = new C(0,0), Pm1 = new C(0,0);
    let D0 = new C(0,0), D1 = new C(0,0), Dm1 = new C(0,0), D2 = new C(0,0), Dm2 = new C(0,0);

    for (let r of resonances){
      const L = Number(r.L ?? r.l ?? r[0]);
      const m = Number(r.m ?? r.M ?? r[1]);
      const massRes = Number(r.mass ?? r.Mass ?? r[2]);
      const widthRes = Number(r.width ?? r.Width ?? r[3]);
      const amp = (typeof r.amplitude !== 'undefined') ? Number(r.amplitude) : ((typeof r.Amplitude !== 'undefined')?Number(r.Amplitude):1.0);
      const x = (typeof xOverride !== 'undefined') ? Number(xOverride) : amp;
      const BW = bwFunc(massRes, widthRes, massKK, x);

      if (L === 0 && m === 0) {
        S0.re += BW.re; S0.im += BW.im;
      } else if (L === 1) {
        if (m === 0) { P0.re += BW.re; P0.im += BW.im; }
        else if (m === 1) { P1.re += BW.re; P1.im += BW.im; }
        else if (m === -1) { Pm1.re += BW.re; Pm1.im += BW.im; }
      } else if (L === 2) {
        if (m === 0) { D0.re += BW.re; D0.im += BW.im; }
        else if (m === 1) { D1.re += BW.re; D1.im += BW.im; }
        else if (m === -1) { Dm1.re += BW.re; Dm1.im += BW.im; }
        else if (m === 2) { D2.re += BW.re; D2.im += BW.im; }
        else if (m === -2) { Dm2.re += BW.re; Dm2.im += BW.im; }
      }
    }

    return { S0, Pm1, P0, P1, Dm2, Dm1, D0, D1, D2 };
  }

  // allMomentsBuilder: compute the 30 H0/H1 moments and return an array
  function allMomentsBuilder(resonances, massKK, xOverride){
    const pw = waveBuilder(resonances, massKK, xOverride);
    const S0 = pw.S0, Pm1 = pw.Pm1, P0 = pw.P0, P1 = pw.P1, Dm2 = pw.Dm2, Dm1 = pw.Dm1, D0 = pw.D0, D1 = pw.D1, D2 = pw.D2;

    const sqrt3  = Math.sqrt(3.0);
    const sqrt5  = Math.sqrt(5.0);
    const sqrt15 = Math.sqrt(15.0);

    const nS0 = norm(S0), nP0 = norm(P0), nP1 = norm(P1), nPm1 = norm(Pm1);
    const nD0 = norm(D0), nD1 = norm(D1), nDm1 = norm(Dm1), nD2 = norm(D2), nDm2 = norm(Dm2);
    function R(a,b){ return realProd(a,b); }

    const H0_00 = 2.0 * ( nS0 + nPm1 + nP0 + nP1 + nDm2 + nDm1 + nD0 + nD1 + nD2 );

    const H0_10 = (4.0 / sqrt3) * R(S0, P0)
                + (8.0 / sqrt15) * R(P0, D0)
                + (4.0 / Math.sqrt(5)) * R(P1, D1)
                + (4.0 / Math.sqrt(5)) * R(Pm1, Dm1);

    const H0_11 = 2.0 * (
                (sqrt3 / 3) * R(S0, P1)
                - (sqrt3 / 3) * R(S0, Pm1)
                + (Math.sqrt(10) / 5) * R(P1, D2)
                - (sqrt15 / 15) * R(P1, D0)
                + (sqrt5 / 5) * R(P0, D1)
                - (sqrt5 / 5) * R(P0, Dm1)
                + (sqrt15 / 15) * R(Pm1, D0)
                - (Math.sqrt(10) / 5) * R(Pm1, Dm2) );

    const H0_20 = (4.0/5.0) * nP0
                - (2.0/5.0) * ( nP1 + nPm1 )
                + (4.0/7.0) * nD0
                + (2.0/7.0) * ( nD1 + nDm1 )
                - (4.0/7.0) * ( nD2 + nDm2 )
                + ((4.0*sqrt5) / 5) * R(S0, D0);

    const H0_21 = 2.0 * (
                (sqrt5 / 5) * R(S0, D1)
                - (sqrt5 / 5) * R(S0, Dm1)
                + (Math.sqrt(3) / 5) * R(P1, P0)
                - (Math.sqrt(3) / 5) * R(P0, Pm1)
                + (Math.sqrt(6) / 7) * R(D2, D1)
                + (1.0 / 7) * R(D1, D0)
                - (1.0 / 7) * R(D0, Dm1)
                - (Math.sqrt(6) / 7) * R(Dm1, Dm2) );

    const H0_22 = 2.0 * ( R( new C(S0.re / sqrt5 - 2.0 * D0.re / 7.0, S0.im / sqrt5 - 2.0 * D0.im / 7.0), new C(D2.re + Dm2.re, D2.im + Dm2.im) ) )
                - 2.0 * (Math.sqrt(6.0)/5.0) * R(P1, Pm1)
                - 2.0 * (Math.sqrt(6.0)/7.0) * R(D1, Dm1);

    const H0_30 = (12.0 * sqrt5 / 35.0) * (
                - R(P1, D1)
                - R(Pm1, Dm1)
                + sqrt3 * R(P0, D0) );

    const H0_31 = 2 * Math.sqrt(5) / 35 * (
                - Math.sqrt(3) * R(P1, D2)
                + 3 * Math.sqrt(2) * R(P1, D0)
                + 2 * Math.sqrt(6) * R(P0, D1)
                - 2 * Math.sqrt(6) * R(P0, Dm1)
                - 3 * Math.sqrt(2) * R(Pm1, D0)
                + Math.sqrt(3) * R(Pm1, Dm2) );

    const H0_32 = (2.0 * Math.sqrt(3.0) / 7.0) * (-Math.sqrt(2.0) * R(P1, Dm1)
                + R(P0, D2) + R(P0, Dm2)
                - Math.sqrt(2.0) * R(Pm1, D1) );

    const H0_33 = (6.0 / 7.0) * ( R(P1, Dm2) - R(Pm1, D2) );

    const H0_40 = (2.0/21.0) * ( nD2 - 4.0*nD1 + 6.0*nD0 - 4.0*nDm1 + nDm2 );

    const H0_41 = (2.0 * sqrt5 / 21.0) * (
                - R(D2, D1)
                + Math.sqrt(6) * R(D1, D0)
                - Math.sqrt(6) * R(D0, Dm1)
                + R(Dm1, Dm2) );

    const H0_42 = ((2 * Math.sqrt(5)) / 21) * ( Math.sqrt(3) * R(D2, D0)
                - 2 * Math.sqrt(2) * R(D1, Dm1) + Math.sqrt(3) * R(D0, Dm2) );

    const H0_43 = (2.0 * Math.sqrt(35) / 21.0) * ( - R(D2, Dm1) + R(D1, Dm2) );

    const H0_44 = ((2*Math.sqrt(70))/21.0) * R(D2, Dm2);

    const H1_00 = 2.0 * ( nS0 + nP0 + nD0
                        - 2.0 * R(P1, Pm1)
                        - 2.0 * R(D1, Dm1)
                        + 2.0 * R(D2, Dm2) );

    const H1_10 = 4 * (
        (Math.sqrt(3) / 3) * R(S0, P0)
        - (Math.sqrt(5) / 5) * R(P1, Dm1)
        + (2 * Math.sqrt(15) / 15) * R(P0, D0)
        - (Math.sqrt(5) / 5) * R(Pm1, D1) );

    const H1_11 = 2 * (
        (Math.sqrt(3) / 3) * R(S0, P1)
        - (Math.sqrt(3) / 3) * R(S0, Pm1)
        - (Math.sqrt(15) / 15) * R(P1, D0)
        + (Math.sqrt(10) / 5) * R(P1, Dm2)
        + (Math.sqrt(5) / 5) * R(P0, D1)
        - (Math.sqrt(5) / 5) * R(P0, Dm1)
        - (Math.sqrt(10) / 5) * R(Pm1, D2)
        + (Math.sqrt(15) / 15) * R(Pm1, D0) );

    const H1_20 = (4.0/5.0) * nP0
                + (4.0/7.0) * R(P1, Pm1)
                + (4.0/7.0) * nD0
                - (4.0/7.0) * R(D1, Dm1)
                - (8.0/7.0) * R(D2, Dm2)
                + (4.0 / Math.sqrt(5.0)) * R(S0, D0);

    const H1_21 = 2 * (
        (Math.sqrt(5) / 5) * R(S0, D1)
        - (Math.sqrt(5) / 5) * R(S0, Dm1)
        + (Math.sqrt(3) / 5) * R(P1, P0)
        - (Math.sqrt(3) / 5) * R(P0, Pm1)
        - (Math.sqrt(6) / 7) * R(D2, Dm1)
        + (1.0 / 7) * R(D1, D0)
        + (Math.sqrt(6) / 7) * R(D1, D2)
        - (1.0 / 7) * R(D0, Dm1) );

    const H1_22 = 2 * (
        (Math.sqrt(5) / 5) * R(S0, D2)
        + (Math.sqrt(5) / 5) * R(S0, Dm2)
        + (Math.sqrt(6) / 10) * nP1
        + (Math.sqrt(6) / 10) * nPm1
        - (2.0 / 7) * R(D2, D0)
        - (2.0 / 7) * R(Dm2, D0)
        + (Math.sqrt(6) / 14) * nD1
        + (Math.sqrt(6) / 14) * nDm1 );

    const H1_30 = (12.0 * sqrt5 / 35.0) * ( R(P1, Dm1)
                + Math.sqrt(3) * R(P0, D0)
                + R(Pm1, D1) );

    const H1_31 = 2 * Math.sqrt(5) / 35 * (
        3 * Math.sqrt(2) * R(P1, D0)
        - Math.sqrt(3) * R(P1, Dm2)
        + 2 * Math.sqrt(6) * R(P0, D1)
        - 2 * Math.sqrt(6) * R(P0, Dm1)
        + Math.sqrt(3) * R(Pm1, D2)
        - 3 * Math.sqrt(2) * R(Pm1, D0) );

    const H1_32 = (2.0 / 7.0) * (
            Math.sqrt(6) * R(P1, D1)
            + Math.sqrt(3) * R(P0, D2)
            + Math.sqrt(3) * R(P0, Dm2)
            + Math.sqrt(6) * R(Pm1, Dm1) );

    const H1_33 = (6.0 / 7.0) * ( R(P1, D2) - R(Pm1, Dm2) );

    const H1_40 = (4.0/21.0) * ( 3.0*nD0 + 4.0*R(D1, Dm1) + R(D2, Dm2) );

    const H1_41 = (2.0*Math.sqrt(5)/21.0) * ( Math.sqrt(6) * R(D0, D1) - Math.sqrt(6) * R(D0, Dm1)
                                        + R(D1, D2) - R(Dm1, Dm2) );

    const H1_42 = (2 * Math.sqrt(5) / 21) * (
        Math.sqrt(3) * R(D2, D0)
        + Math.sqrt(2) * nDm1
        + Math.sqrt(3) * R(D0, Dm2)
    );

    const H1_43 = ((2.0*Math.sqrt(35))/21.0) * ( R(D2, D1) - R(Dm1, Dm2) );

    const H1_44 = (Math.sqrt(70)/21) * ( nD2 + nDm2 );

    return [H0_00, H0_10, H0_11, H0_20, H0_21, H0_22, H0_30, H0_31, H0_32, H0_33,
            H0_40, H0_41, H0_42, H0_43, H0_44,
            H1_00, H1_10, H1_11, H1_20, H1_21, H1_22, H1_30, H1_31, H1_32, H1_33,
            H1_40, H1_41, H1_42, H1_43, H1_44];
  }

  // export
  window.bwFunc = bwFunc;
  window.waveBuilder = waveBuilder;
  window.allMomentsBuilder = allMomentsBuilder;

})();
