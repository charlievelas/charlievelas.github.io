// app.js
// Client-side plotting for EtaPi moments.
// - Generates dynamic plots for all moments H^0(LM) and -H^1(LM) for L=0..4, M=0..L
// - Generates BA plots (BA 4pi, BA along y, binned BA).
// - Data source:
//     * If a client-side compute function `computeMomentsClient(params)` is provided (pure JS port or WASM wrapper),
//       the page will call it and use its results (preferred).
//     * Otherwise the demo will fetch existing data files moment0.txt, moment1.txt, BA.txt from the repo's
//       EtaPi_moments folder and plot them. This keeps the plotting dynamic (no PNGs).
//
// Next step: replace the "fallback fetch" by a true JS/WASM compute routine so the outputs depend on the inputs.

(function () {
  const startBtn = document.getElementById('startBtn');
  const status = document.getElementById('status');
  const downloads = document.getElementById('downloads');
  const plotsContainer = document.getElementById('plots');

  // Column ordering in the moment files (per PHP/gplot):
  // columns: m, 00, 10,11, 20,21,22, 30,31,32,33, 40,41,42,43,44
  const LM_ORDER = [
    {L:0,M:0},                 // col 1 (index 1 in file, but index 0 here after m)
    {L:1,M:0}, {L:1,M:1},
    {L:2,M:0}, {L:2,M:1}, {L:2,M:2},
    {L:3,M:0}, {L:3,M:1}, {L:3,M:2}, {L:3,M:3},
    {L:4,M:0}, {L:4,M:1}, {L:4,M:2}, {L:4,M:3}, {L:4,M:4}
  ];

  // Raw URLs (adjust branch if needed)
  const rawBase = 'https://raw.githubusercontent.com/charlievelas/charlievelas.github.io/main/EtaPi_moments/';
  const urls = {
    moment0: rawBase + 'moment0.txt',
    moment1: rawBase + 'moment1.txt',
    ba:      rawBase + 'BA.txt'
  };

  // Build plot boxes for each LM in the same sequence as gnuplot.txt
  function buildPlotBoxes() {
    plotsContainer.innerHTML = '';
    // H moments
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const {L, M} = LM_ORDER[idx];
      const id = `pH${L}${M}`;
      const box = document.createElement('div');
      box.className = 'plotBox';
      box.id = id;
      const title = document.createElement('div');
      title.innerHTML = `<strong>H^0(${L}${M}) & -H^1(${L}${M})</strong>`;
      title.style.marginBottom = '6px';
      box.appendChild(title);
      plotsContainer.appendChild(box);
    }
    // BA plots (two)
    const baBox1 = document.createElement('div');
    baBox1.className = 'plotBox';
    baBox1.id = 'pBA';
    baBox1.innerHTML = '<strong>Beam Asymmetry Integrated (BA 4π and BA_y)</strong>';
    plotsContainer.appendChild(baBox1);

    const baBox2 = document.createElement('div');
    baBox2.className = 'plotBox';
    baBox2.id = 'pBAy';
    baBox2.innerHTML = '<strong>Beam Asymmetry Binned / Additional</strong>';
    plotsContainer.appendChild(baBox2);
  }

  // Parse whitespace-delimited numeric table into array of rows (numbers)
  function parseTable(text) {
    const rows = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length && !l.startsWith('#'));
    return rows.map(line => line.split(/\s+/).map(Number));
  }

  // Plot all H^0(LM) and -H^1(LM)
  function plotAllMoments(moment0Rows, moment1Rows) {
    // moment0Rows/moment1Rows: arrays of number arrays
    // both expected to have same m column length; if not, use intersection by m or fallback to moment0 m
    const m0 = moment0Rows.map(r => r[0]);

    // For each LM in LM_ORDER, pick corresponding column index in file:
    // File indices (one-based): 1=m, 2=00, 3=10,4=11,5=20 ... so mapping colIndex = index_in_LM_ORDER + 1
    for (let idx = 0; idx < LM_ORDER.length; idx++) {
      const colIndex = idx + 1; // zero-based row array: 0=m, 1=00, 2=10...
      const L = LM_ORDER[idx].L;
      const M = LM_ORDER[idx].M;
      const divId = `pH${L}${M}`;
      const y0 = moment0Rows.map(r => r[colIndex]);
      const y1 = moment1Rows.map(r => r[colIndex]);
      const trace0 = { x: m0, y: y0, mode: 'lines', name: `H^0(${L}${M})`, line: { color: 'red' } };
      const trace1 = { x: m0, y: y1.map(v => -v), mode: 'lines', name: `-H^1(${L}${M})`, line: { color: 'blue' } };

      const layout = { margin: { t: 24, b: 36, l: 44, r: 8 }, xaxis: { title: 'm (GeV)' }, yaxis: { title: `H(${L}${M})` } };
      Plotly.newPlot(divId, [trace0, trace1], layout, { responsive: true });
    }
  }

  // Plot BA outputs (BA.txt expected format: m ba4pi bay ba_binned)
  function plotBA(baRows) {
    const m = baRows.map(r => r[0]);
    const ba4pi = baRows.map(r => r[1]);
    const bay = baRows.map(r => r[2]);
    const ba_binned = baRows.map(r => (r.length > 3 ? r[3] : NaN));

    Plotly.newPlot('pBA', [
      { x: m, y: ba4pi, mode: 'lines', name: 'BA 4π', line: { color: 'red' } },
      { x: m, y: bay, mode: 'lines', name: 'BA along y', line: { color: 'blue' } }
    ], { margin: { t: 24 }, xaxis: { title: 'm (GeV)' } });

    Plotly.newPlot('pBAy', [
      { x: m, y: ba_binned, mode: 'lines', name: 'BA binned', line: { color: 'green' } }
    ], { margin: { t: 24 }, xaxis: { title: 'm (GeV)' } });
  }

  // Create download links for produced data
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

  // Fallback: fetch precomputed numeric files from repo and plot them (no PNGs)
  async function fetchAndPlotFallback() {
    status.textContent = 'Fetching precomputed numeric files...';
    try {
      const [r0, r1, rb] = await Promise.all([fetch(urls.moment0), fetch(urls.moment1), fetch(urls.ba)]);
      if (!r0.ok || !r1.ok || !rb.ok) {
        console.warn('One or more fetches failed', r0.status, r1.status, rb.status);
        status.textContent = 'Error fetching numeric files (see console).';
        return;
      }
      const [txt0, txt1, txtb] = await Promise.all([r0.text(), r1.text(), rb.text()]);
      const rows0 = parseTable(txt0);
      const rows1 = parseTable(txt1);
      const rowsb = parseTable(txtb);

      plotAllMoments(rows0, rows1);
      plotBA(rowsb);

      downloads.innerHTML = '';
      downloads.appendChild(makeDownloadLink('moment0.txt', txt0));
      downloads.appendChild(makeDownloadLink('moment1.txt', txt1));
      downloads.appendChild(makeDownloadLink('BA.txt', txtb));

      status.textContent = 'Plotted precomputed numeric files.';
    } catch (err) {
      console.error(err);
      status.textContent = 'Error fetching/plotting: ' + err.message;
    }
  }

  // Main entry: compute & plot.
  // If a real client-side compute routine exists, call it.
  // Expect computeMomentsClient(params) -> { moment0Text, moment1Text, baText }
  async function computeAndPlot() {
    // Read form parameters
    const params = {
      Eg: Number(document.getElementById('Eg').value),
      t:  Number(document.getElementById('t').value),
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

    status.textContent = 'Starting computation...';

    // If user (or a follow-up task) supplies a client-side compute function, call it.
    // The function should return an object with { moment0Text, moment1Text, baText } (strings like files),
    // or arrays (moment0Rows, moment1Rows, baRows).
    if (typeof window.computeMomentsClient === 'function') {
      status.textContent = 'Running client-side compute routine...';
      try {
        const out = await window.computeMomentsClient(params);
        // out may return parsed rows or raw text; handle both
        let rows0, rows1, rowsb;
        if (out.moment0Text) rows0 = parseTable(out.moment0Text);
        else if (out.moment0Rows) rows0 = out.moment0Rows;
        if (out.moment1Text) rows1 = parseTable(out.moment1Text);
        else if (out.moment1Rows) rows1 = out.moment1Rows;
        if (out.baText) rowsb = parseTable(out.baText);
        else if (out.baRows) rowsb = out.baRows;

        if (!rows0 || !rows1) {
          status.textContent = 'Computation returned incomplete data; falling back to precomputed files.';
          await fetchAndPlotFallback();
          return;
        }

        plotAllMoments(rows0, rows1);
        if (rowsb) plotBA(rowsb);

        downloads.innerHTML = '';
        if (out.moment0Text) downloads.appendChild(makeDownloadLink('moment0.txt', out.moment0Text));
        if (out.moment1Text) downloads.appendChild(makeDownloadLink('moment1.txt', out.moment1Text));
        if (out.baText) downloads.appendChild(makeDownloadLink('BA.txt', out.baText));

        status.textContent = 'Plotted results from client-side compute routine.';
        return;
      } catch (err) {
        console.error('Client compute failed:', err);
        status.textContent = 'Client compute failed (see console). Falling back to precomputed files.';
        await fetchAndPlotFallback();
        return;
      }
    }

    // Default: fallback to fetching the repo files and plotting them dynamically
    status.textContent = 'No client compute routine found; plotting precomputed numeric files.';
    await fetchAndPlotFallback();
  }

  // Initialize UI
  buildPlotBoxes();
  status.textContent = 'Ready.';

  startBtn.addEventListener('click', async function () {
    status.textContent = 'Preparing plots...';
    // Rebuild plot boxes to ensure a fresh set
    buildPlotBoxes();
    await computeAndPlot();
  });

})();
