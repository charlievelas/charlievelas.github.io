// Improved app.js — better fetch logic and clearer error reporting.
// Tries relative paths first (works when hosted on GitHub Pages) and falls back to raw.githubusercontent.
// Also logs HTTP status codes and displays helpful status messages.

(function () {
  const startBtn = document.getElementById('startBtn');
  const status = document.getElementById('status');
  const downloads = document.getElementById('downloads');
  const plotsContainer = document.getElementById('plots');

  const LM_ORDER = [
    {L:0,M:0},
    {L:1,M:0}, {L:1,M:1},
    {L:2,M:0}, {L:2,M:1}, {L:2,M:2},
    {L:3,M:0}, {L:3,M:1}, {L:3,M:2}, {L:3,M:3},
    {L:4,M:0}, {L:4,M:1}, {L:4,M:2}, {L:4,M:3}, {L:4,M:4}
  ];

  // Prefer relative repo paths (will work with GitHub Pages). If those fail, fall back to raw.githubusercontent
  const relativeBase = './EtaPi_moments/';
  const rawBase = 'https://raw.githubusercontent.com/charlievelas/charlievelas.github.io/main/EtaPi_moments/';

  const fileNames = {
    moment0: 'moment0.txt',
    moment1: 'moment1.txt',
    ba:      'BA.txt'
  };

  function buildPlotBoxes() {
    plotsContainer.innerHTML = '';
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
    const baBox1 = document.createElement('div'); baBox1.className='plotBox'; baBox1.id='pBA';
    baBox1.innerHTML = '<strong>Beam Asymmetry Integrated (BA 4π and BA_y)</strong>';
    plotsContainer.appendChild(baBox1);
    const baBox2 = document.createElement('div'); baBox2.className='plotBox'; baBox2.id='pBAy';
    baBox2.innerHTML = '<strong>Beam Asymmetry Binned / Additional</strong>';
    plotsContainer.appendChild(baBox2);
  }

  function parseTable(text) {
    const rows = text.split(/\r?\n/).map(l => l.trim()).filter(l => l.length && !l.startsWith('#'));
    return rows.map(line => line.split(/\s+/).map(Number));
  }

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

  // Attempt to fetch a file first via relative path, then raw.githubusercontent fallback.
  async function fetchWithFallback(fileName) {
    const relUrl = relativeBase + fileName;
    const rawUrl = rawBase + fileName;
    // Try relative
    try {
      const r = await fetch(relUrl, {cache:'no-cache'});
      if (r.ok) {
        console.log(`Fetched ${fileName} via relative path (${relUrl}) status=${r.status}`);
        return { ok:true, text: await r.text(), url: relUrl, status: r.status };
      } else {
        console.warn(`Relative fetch failed for ${fileName}: ${relUrl} status=${r.status}`);
      }
    } catch (err) {
      console.warn(`Relative fetch error for ${fileName}:`, err.message || err);
    }
    // Try raw.githubusercontent fallback
    try {
      const r2 = await fetch(rawUrl, {cache:'no-cache'});
      if (r2.ok) {
        console.log(`Fetched ${fileName} via raw.githubusercontent (${rawUrl}) status=${r2.status}`);
        return { ok:true, text: await r2.text(), url: rawUrl, status: r2.status };
      } else {
        console.warn(`Raw fetch failed for ${fileName}: ${rawUrl} status=${r2.status}`);
        return { ok:false, status: r2.status, url: rawUrl };
      }
    } catch (err) {
      console.error(`Raw fetch error for ${fileName}:`, err.message || err);
      return { ok:false, error: err, url: rawUrl };
    }
  }

  async function fetchAndPlotFallback() {
    status.textContent = 'Fetching numeric files (relative → raw fallback)...';
    try {
      const res0 = await fetchWithFallback(fileNames.moment0);
      const res1 = await fetchWithFallback(fileNames.moment1);
      const resb = await fetchWithFallback(fileNames.ba);

      const failed = [];
      if (!res0.ok) failed.push({name:fileNames.moment0, url:res0.url, status:res0.status});
      if (!res1.ok) failed.push({name:fileNames.moment1, url:res1.url, status:res1.status});
      if (!resb.ok) failed.push({name:fileNames.ba, url:resb.url, status:resb.status});

      if (failed.length) {
        console.warn('Fetch failures:', failed);
        status.textContent = 'Error fetching files. Check console/Network for details. Missing: ' + failed.map(f=>f.name+'('+f.status+')').join(', ');
        // show exact URLs in console
        failed.forEach(f => console.warn('Failed file:', f));
        return;
      }

      const rows0 = parseTable(res0.text);
      const rows1 = parseTable(res1.text);
      const rowsb = parseTable(resb.text);

      plotAllMoments(rows0, rows1);
      plotBA(rowsb);

      downloads.innerHTML = '';
      downloads.appendChild(makeDownloadLink('moment0.txt', res0.text));
      downloads.appendChild(makeDownloadLink('moment1.txt', res1.text));
      downloads.appendChild(makeDownloadLink('BA.txt', resb.text));
      status.textContent = `Plotted files from ${res0.url}, ${res1.url}, ${resb.url}`;
    } catch (err) {
      console.error('fetchAndPlotFallback error:', err);
      status.textContent = 'Error fetching/plotting: ' + (err.message || err);
    }
  }

  async function computeAndPlot() {
    // read minimal form parameters (not used by fallback)
    // If you implement window.computeMomentsClient(params) it will be used instead.
    const params = { Eg: Number(document.getElementById('Eg').value), t: Number(document.getElementById('t').value), dm: Number(document.getElementById('dm').value) };
    status.textContent = 'Starting compute & plot...';

    if (typeof window.computeMomentsClient === 'function') {
      status.textContent = 'Calling client-side compute routine...';
      try {
        const out = await window.computeMomentsClient(params);
        let rows0, rows1, rowsb;
        if (out.moment0Text) rows0 = parseTable(out.moment0Text);
        else if (out.moment0Rows) rows0 = out.moment0Rows;
        if (out.moment1Text) rows1 = parseTable(out.moment1Text);
        else if (out.moment1Rows) rows1 = out.moment1Rows;
        if (out.baText) rowsb = parseTable(out.baText);
        else if (out.baRows) rowsb = out.baRows;

        if (!rows0 || !rows1) {
          status.textContent = 'Client compute returned incomplete data; falling back to file fetch.';
          await fetchAndPlotFallback();
          return;
        }
        plotAllMoments(rows0, rows1);
        if (rowsb) plotBA(rowsb);
        downloads.innerHTML = '';
        if (out.moment0Text) downloads.appendChild(makeDownloadLink('moment0.txt', out.moment0Text));
        if (out.moment1Text) downloads.appendChild(makeDownloadLink('moment1.txt', out.moment1Text));
        if (out.baText) downloads.appendChild(makeDownloadLink('BA.txt', out.baText));
        status.textContent = 'Plotted results from client compute routine.';
        return;
      } catch (err) {
        console.error('client compute error:', err);
        status.textContent = 'Client compute failed. Falling back to file fetch.';
      }
    }

    // fallback
    await fetchAndPlotFallback();
  }

  buildPlotBoxes();
  status.textContent = 'Ready.';
  startBtn.addEventListener('click', async function () {
    buildPlotBoxes();
    await computeAndPlot();
  });

})();
