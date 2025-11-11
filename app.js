// app.js
// Client-side plotting for EtaPi moments — no external files required.
// This version ALWAYS generates numeric data client-side from the form inputs
// (a simple physics-inspired synthetic generator). Later you can replace
// computeMomentsClient() with a full JS port or a WASM wrapper that returns
// the same output format.
//
// Output format produced by computeMomentsClient:
//  - moment0Text: text of moment0.txt (m + 15 columns: 00,10,11,20,21,22,30,...,44)
//  - moment1Text: text of moment1.txt (same columns, sign convention per original code)
//  - baText:      text of BA.txt (m, ba4pi, bay, ba_binned)
//
// The plotting code expects row-based numeric arrays or the text and will plot them.

(function () {
  const startBtn = document.getElementById('startBtn');
  const statusEl = document.getElementById('status');
  const downloads = document.getElementById('downloads');
  const plotsContainer = document.getElementById('plots');

  // LM order used by gnuplot.txt / original PHP: columns after m
  const LM_ORDER = [
    {L:0,M:0},                 // 00
    {L:1,M:0}, {L:1,M:1},      // 10,11
    {L:2,M:0}, {L:2,M:1}, {L:2,M:2}, // 20,21,22
    {L:3,M:0}, {L:3,M:1}, {L:3,M:2}, {L:3,M:3}, // 30..33
    {L:4,M:0}, {L:4,M:1}, {L:4,M:2}, {L:4,M:3}, {L:4,M:4} // 40..44
  ];

  function buildPlotBoxes() {
    plotsContainer.innerHTML = '';
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const {L, M} = LM_ORDER[idx];
      const id = `pH${L}${M}`;
      const box = document.createElement('div');
      box.className = 'plotBox';
      box.id = id;
      const title = document.createElement('div');
      title.innerHTML = `<strong>H^0(${L}${M}) &nbsp;&nbsp; &nbsp;-H^1(${L}${M})</strong>`;
      title.style.marginBottom = '6px';
      box.appendChild(title);
      plotsContainer.appendChild(box);
    }
    const baBox1 = document.createElement('div'); baBox1.className='plotBox'; baBox1.id='pBA';
    baBox1.innerHTML = '<strong>Beam Asymmetry Integrated (BA 4π and BA_y)</strong>';
    plotsContainer.appendChild(baBox1);
    const baBox2 = document.createElement('div'); baBox2.className='plotBox'; baBox2.id='pBAy';
    baBox2.innerHTML = '<strong>Beam Asymmetry Binned / Additional</strong>';
    plotsContainer.appendChild(baBox2);
  }

  // parse text table to rows of numbers
  function parseTable(text) {
    const rows = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length && !l.startsWith('#'));
    return rows.map(line => line.split(/\s+/).map(Number));
  }

  // plotting helpers (Plotly)
  function plotAllMoments(moment0Rows, moment1Rows) {
    const m0 = moment0Rows.map(r => r[0]);
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const colIndex = idx + 1; // 0=m, 1=00, 2=10...
      const {L, M} = LM_ORDER[idx];
      const divId = `pH${L}${M}`;
      const y0 = moment0Rows.map(r => r[colIndex]);
      const y1 = moment1Rows.map(r => r[colIndex]);
      const trace0 = { x: m0, y: y0, mode: 'lines', name: `H^0(${L}${M})`, line: { color: 'red' } };
      const trace1 = { x: m0, y: y1.map(v => -v), mode: 'lines', name: `-H^1(${L}${M})`, line: { color: 'blue' } };
      const layout = { margin: { t: 24, b: 36, l: 44, r: 8 }, xaxis: { title: 'm (GeV)' }, yaxis: { title: `H(${L}${M})` } };
      Plotly.newPlot(divId, [trace0, trace1], layout, { responsive: true });
    }
  }

  function plotBA(baRows) {
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

  // --------- Synthetic generator (replaceable) ----------
  // This function creates moment0/moment1/ba data from input params.
  // It is a simple physics-inspired generator using resonances (Breit-Wigner magnitudes)
  // and shape factors for different L,M to produce plausible-looking curves.
  //
  // Signature: async function computeMomentsClient(params) -> { moment0Text, moment1Text, baText, moment0Rows, moment1Rows, baRows }
  async function computeMomentsClient(params) {
    // constants (same as main.h)
    const MP = 0.93827203, MPI = 0.13957061, META = 0.547682;
    // collect resonance params from form (m,g). Use defaults if missing.
    const res = [
      { m: params.m0 || 0.98, g: params.g0 || 0.075, N: 1.0, del: 1.0 },
      { m: params.m1 || 1.564, g: params.g1 || 0.492, N: 1.0, del: -5.0 },
      { m: params.m2 || 1.318, g: params.g2 || 0.107, N: 1.0, del: -2.0 },
      { m: params.m2p || 1.722, g: params.g2p || 0.247, N: 1.0, del: -2.0 }
    ];

    // build mass grid
    let dm = Number(params.dm) || 0.02;
    if (dm <= 0 || !isFinite(dm)) dm = 0.02;
    const mMin = META + MPI;
    const mMax = 2.0;
    const masses = [];
    for (let m = mMin; m <= mMax + 1e-12; m += dm) masses.push(parseFloat(m.toFixed(6)));

    // helper: Breit-Wigner magnitude (real positive) using simple form: |BW| = (m*Gamma) / sqrt((mR^2 - m^2)^2 + (mR*Gamma)^2)
    function bwMag(mass, mR, gR) {
      const num = mR * gR;
      const den = Math.sqrt(Math.pow(mR*mR - mass*mass, 2) + Math.pow(mR*gR, 2));
      return num / Math.max(den, 1e-12);
    }

    // energy-dependent scale to mimic s^alpha(t) behaviour (rough)
    const Eg = Number(params.Eg) || 8.5;
    const tval = Number(params.t) || -0.1;
    const sScale = Math.pow(1 + (Eg/10), 0.5) * Math.exp(-Math.abs(tval) * 0.5);

    // beam polarization factor proxy (for H^1 relative to H^0)
    const Pg = 1.0; // assume P_gamma = 1 for scaling H1 as original used -2/Pg etc
    const asymFactor = 0.25 + 0.25*Math.tanh((Eg-6.0)/3.0) * Math.exp(-Math.abs(tval)*0.8);

    // compute H^0(LM) and H^1(LM) arrays
    const moment0Rows = []; // rows: [m, col1..col15]
    const moment1Rows = [];
    const baRows = []; // rows: [m, ba4pi, bay, ba_binned]

    for (const m of masses) {
      // build resonance-weighted amplitude sum and a shape multiplier that depends on L and M
      let ampSum = 0;
      const bwVals = res.map(r => bwMag(m, r.m, r.g));
      for (let i=0; i<res.length; i++) ampSum += bwVals[i] * res[i].N;

      // baseline intensity for H^0(00)
      const H00 = Math.max(0, sScale * ampSum * (1 + 0.6*Math.exp(-Math.pow(m - 1.0, 2)/0.08)));

      // produce columns for LM_ORDER using shape factors
      const cols0 = [];
      const cols1 = [];
      for (let idx = 0; idx < LM_ORDER.length; idx++) {
        const {L,M} = LM_ORDER[idx];
        // shape: higher L suppressed; M modulates oscillation
        const Lscale = Math.exp(-0.8 * L);
        const Mscale = 1 + 0.25 * M;
        // add small oscillatory behavior to mimic interference
        const osc = 1 + 0.15 * Math.sin(3.0 * m * (1 + 0.15*L + 0.05*M));
        // column amplitude
        const val0 = H00 * Lscale * Mscale * osc * (0.6 + 0.4 * Math.cos((L+1)*m*1.2));
        // H1 is proportional to H0 but with sign/phase changes; include asymFactor and an L-dependent sign
        const sign = ((L+M) % 2 === 0) ? 1 : -1;
        const val1 = sign * val0 * asymFactor * (0.8 + 0.2*Math.sin(m*2.5));
        cols0.push(val0);
        cols1.push(val1);
      }

      // Prepare BA values: ba4pi = -H1(00)/H0(00) (mimic original code: -creal(mom1/mom0))
      const ba4pi = (cols0[0] !== 0) ? -cols1[0]/cols0[0] : 0;
      // bay: a different asymmetry (use a local computed version)
      const bay = Math.tanh(0.8 * (Math.sin(2.0*m) * 0.2 + ba4pi * 0.7));
      // ba_binned: simulate binned BA (smoothed)
      const ba_binned = Math.tanh(0.9 * (ba4pi * 0.8 + 0.1*Math.cos(m*5)));

      moment0Rows.push([m, ...cols0]);
      moment1Rows.push([m, ...cols1]);
      baRows.push([m, ba4pi, bay, ba_binned]);
    }

    // convert to text format like the original files (space separated)
    function rowsToText(rows) {
      return rows.map(r => r.map(v => Number.isFinite(v) ? v.toPrecision(8) : 'NaN').join(' ')).join('\n');
    }

    const moment0Text = rowsToText(moment0Rows);
    const moment1Text = rowsToText(moment1Rows);
    const baText = rowsToText(baRows);

    return {
      moment0Text, moment1Text, baText,
      moment0Rows, moment1Rows, baRows
    };
  }
  // --------- end synthetic generator ----------

  // Main compute & plot flow
  async function computeAndPlot() {
    statusEl.textContent = 'Computing moments in browser...';
    // read params from form
    const params = {
      Eg: Number(document.getElementById('Eg').value),
      t: Number(document.getElementById('t').value),
      dm: Number(document.getElementById('dm').value),
      FR: Array.from(document.getElementsByName('FR')).find(r => r.checked).value,
      m0: Number(document.getElementById('m0').value),
      g0: Number(document.getElementById('g0').value),
      m1: Number(document.getElementById('m1').value),
      g1: Number(document.getElementById('g1').value),
      m2: Number(document.getElementById('m2').value),
      g2: Number(document.getElementById('g2').value),
      m2p: Number(document.getElementById('m2p').value),
      g2p: Number(document.getElementById('g2p').value)
    };

    try {
      // If you later replace computeMomentsClient with a real port/WASM, it will be used here.
      let out;
      if (typeof window.computeMomentsClient === 'function' && window.computeMomentsClient !== computeMomentsClient) {
        // call external provided compute routine
        statusEl.textContent = 'Running provided compute routine...';
        out = await window.computeMomentsClient(params);
      } else {
        // use built-in synthetic generator
        out = await computeMomentsClient(params);
      }

      // parse and plot
      const rows0 = out.moment0Rows || parseTable(out.moment0Text || '');
      const rows1 = out.moment1Rows || parseTable(out.moment1Text || '');
      const rowsb = out.baRows || parseTable(out.baText || '');

      plotAllMoments(rows0, rows1);
      if (rowsb && rowsb.length) plotBA(rowsb);

      // provide downloads
      downloads.innerHTML = '';
      downloads.appendChild(makeDownloadLink('moment0.txt', out.moment0Text));
      downloads.appendChild(makeDownloadLink('moment1.txt', out.moment1Text));
      downloads.appendChild(makeDownloadLink('BA.txt', out.baText));

      statusEl.textContent = 'Done — plots generated in browser (no external files required).';
    } catch (err) {
      console.error('Computation/plotting error', err);
      statusEl.textContent = 'Error during computation: ' + (err && err.message ? err.message : String(err));
    }
  }

  // initialize
  buildPlotBoxes();
  statusEl.textContent = 'Ready. Click Start to compute and plot.';
  startBtn.addEventListener('click', async () => {
    buildPlotBoxes(); // clear previous plots
    await computeAndPlot();
  });

})();
