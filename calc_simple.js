// calc.js — simple client-side computation for the demo page.
// Computes square and square root of the user-supplied number and
// offers a downloadable text file with the results.

(function () {
  const elNum = document.getElementById('num');
  const btnCompute = document.getElementById('computeBtn');
  const btnReset = document.getElementById('resetBtn');
  const resultsDiv = document.getElementById('results');
  const downloadLink = document.getElementById('downloadLink');

  function formatNumber(n) {
    // Show up to 10 decimal places but trim trailing zeros
    if (!isFinite(n)) return String(n);
    return Number.parseFloat(n).toPrecision(12).replace(/(?:\.0+|(\.\d+?)0+)$/, '$1');
  }

  function computeAndDisplay() {
    const raw = elNum.value;
    if (raw === '') {
      resultsDiv.innerHTML = '<div class="error">Please provide a number.</div>';
      downloadLink.style.display = 'none';
      return;
    }

    const x = Number(raw);
    if (Number.isNaN(x)) {
      resultsDiv.innerHTML = '<div class="error">Invalid number.</div>';
      downloadLink.style.display = 'none';
      return;
    }

    const square = x * x;
    let sqrt;
    let sqrtMessage;

    if (x < 0) {
      // Square root of negative real number is not a real number
      sqrt = NaN;
      sqrtMessage = 'Square root (real): not defined for negative input';
    } else {
      sqrt = Math.sqrt(x);
      sqrtMessage = formatNumber(sqrt);
    }

    const html = [
      `<strong>Input:</strong> ${formatNumber(x)}<br>`,
      `<strong>Square:</strong> ${formatNumber(square)}<br>`,
      `<strong>Square root:</strong> ${sqrtMessage}`
    ].join('\n');

    resultsDiv.innerHTML = html;

    // Prepare text to download
    const text = [
      `input: ${x}`,
      `square: ${square}`,
      `sqrt: ${Number.isNaN(sqrt) ? 'NaN' : sqrt}`
    ].join('\n');

    const blob = new Blob([text], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    downloadLink.href = url;
    downloadLink.style.display = 'inline';
    downloadLink.textContent = 'Download results';
  }

  btnCompute.addEventListener('click', computeAndDisplay);

  btnReset.addEventListener('click', function () {
    elNum.value = '4';
    resultsDiv.textContent = 'Results will appear here after you click Compute.';
    downloadLink.style.display = 'none';
  });

  // Allow pressing Enter while focused on the number field to compute
  elNum.addEventListener('keydown', function (ev) {
    if (ev.key === 'Enter') {
      ev.preventDefault();
      computeAndDisplay();
    }
  });

  // Initialize (optional)
  resultsDiv.textContent = 'Results will appear here after you click Compute.';
})();
