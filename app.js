// app.js
// Basic client-side UI + plotting for EtaPi moments demo.
// For this quick test we fetch precomputed output files from the repo
// and plot them using Plotly. Next step: replace fetch with a JS port
// of the computation so results are produced in-browser.

(function () {
  const startBtn = document.getElementById('startBtn');
  const status = document.getElementById('status');
  const downloadsDiv = document.getElementById('downloads');

  // Raw file URLs (public GitHub raw content). Adjust branch/path if needed.
  const rawBase = 'https://raw.githubusercontent.com/charlievelas/charlievelas.github.io/main/EtaPi_moments/';
  const urls = {
    moment0: rawBase + 'moment0.txt',
    moment1: rawBase + 'moment1.txt',
    ba:      rawBase + 'BA.txt'
  };

  function parseTable(text) {
    // Returns array of rows: each row is array of numbers.
    const rows = text.split(/\r?\n/).map(line => line.trim()).filter(line => line && !line.startsWith('#'));
    return rows.map(line => line.split(/\s+/).map(Number));
  }

  function plotMomentPair(moment0, moment1) {
    // expect each: rows with columns: m, H00, H10, H11, ...
    const m0 = moment0.map(r => r[0]);
    const H00_0 = moment0.map(r => r[1]);
    const H10_0 = moment0.map(r => r[2]);
    const H11_0 = moment0.map(r => r[3]);

    const m1 = moment1.map(r => r[0]);
    const H00_1 = moment1.map(r => r[1]);
    const H10_1 = moment1.map(r => r[2]);

    // Plot H00 (moment0 vs moment1 overlay)
    Plotly.newPlot('plotH00', [
      { x: m0, y: H00_0, mode: 'lines', name: 'H^0(00)', line:{color:'red'} },
      { x: m1, y: H00_1.map(v=> -v), mode: 'lines', name: '-H^1(00)', line:{color:'blue'} }
    ], { title: 'H^0(00) and -H^1(00)', xaxis:{title:'m (GeV)'} });

    // Plot H10
    Plotly.newPlot('plotH10', [
      { x: m0, y: H10_0, mode: 'lines', name: 'H^0(10)', line:{color:'red'} },
      { x: m1, y: H10_1.map(v=> -v), mode: 'lines', name: '-H^1(10)', line:{color:'blue'} }
    ], { title: 'H^0(10) and -H^1(10)', xaxis:{title:'m (GeV)'} });
  }

  function plotBA(baRows) {
    // BA.txt format (per PHP): columns: m, ba4pi, bay, ba_binned
    const m = baRows.map(r => r[0]);
    const ba4pi = baRows.map(r => r[1]);
    const bay = baRows.map(r => r[2]);
    const ba_binned = baRows.map(r => r[3] || NaN);

    Plotly.newPlot('plotBA', [
      { x: m, y: ba4pi, mode:'lines', name:'BA 4π', line:{color:'red'} },
      { x: m, y: bay, mode:'lines', name:'BA along y', line:{color:'blue'} },
      { x: m, y: ba_binned, mode:'lines', name:'BA binned', line:{color:'green'} }
    ], { title: 'Beam Asymmetries', xaxis:{title:'m (GeV)'} });
  }

  function makeDownloadLink(name, text) {
    const blob = new Blob([text], {type:'text/plain'});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = name;
    a.textContent = 'Download ' + name;
    a.style.marginRight = '12px';
    return a;
  }

  async function loadAndPlot() {
    try {
      status.textContent = 'Fetching moment0.txt ...';
      const [m0res, m1res, bares] = await Promise.all([fetch(urls.moment0), fetch(urls.moment1), fetch(urls.ba)]);
      if (!m0res.ok || !m1res.ok || !bares.ok) {
        status.textContent = 'Error fetching files. Check file URLs or network (look in console).';
        console.error('Fetch statuses', m0res.status, m1res.status, bares.status);
        return;
      }
      status.textContent = 'Parsing files ...';
      const [m0txt, m1txt, batxt] = await Promise.all([m0res.text(), m1res.text(), bares.text()]);
      const moment0 = parseTable(m0txt);
      const moment1 = parseTable(m1txt);
      const ba = parseTable(batxt);

      status.textContent = 'Plotting ...';
      plotMomentPair(moment0, moment1);
      plotBA(ba);

      downloadsDiv.innerHTML = '';
      downloadsDiv.appendChild(makeDownloadLink('moment0.txt', m0txt));
      downloadsDiv.appendChild(makeDownloadLink('moment1.txt', m1txt));
      downloadsDiv.appendChild(makeDownloadLink('BA.txt', batxt));

      status.textContent = 'Done: plotted fetched files. To compute results client-side, ask me to port the C code into JS (next step).';
    } catch (err) {
      console.error(err);
      status.textContent = 'Error: ' + err.message;
    }
  }

  startBtn.addEventListener('click', function () {
    status.textContent = 'Starting fetch/plot...';
    loadAndPlot();
  });

})();
