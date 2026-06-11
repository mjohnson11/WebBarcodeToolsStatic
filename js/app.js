// ── Constants ─────────────────────────────────────────────────────────────────

const BARSEQ_DEFAULT = 'GATGTCCACGAGGTCTCTNNNNNNNNNNNNNNNNNNNNCGTACGCTGCAGGTCGAC';
const BC_COLORS      = ['#27ae60', '#2980b9', '#d35400', '#8e44ad'];

// ── State ─────────────────────────────────────────────────────────────────────

let workerInstance = null;
let currentResult  = null;  // {sortedCounts, readSamples, failSamples, stats, regions}
let constructSeq   = BARSEQ_DEFAULT;

// ── DOM refs ──────────────────────────────────────────────────────────────────

const constructInput    = document.getElementById('construct-input');
const constructTextWrap = document.getElementById('construct-text-wrap');
const constructFastaWrap= document.getElementById('construct-fasta-wrap');
const fastaUpload       = document.getElementById('fasta-upload');
const fastaHint         = document.getElementById('fasta-hint');
const fastqUpload       = document.getElementById('fastq-upload');
const parseBtn          = document.getElementById('parse-btn');
const parseHint         = document.getElementById('parse-hint');
const progressWrap      = document.getElementById('progress-wrap');
const progressFill      = document.getElementById('progress-fill');
const progressLabel     = document.getElementById('progress-label');
const resultsEl         = document.getElementById('results');
const statCards         = document.getElementById('stat-cards');
const regionStatsSection= document.getElementById('region-stats-section');
const regionStatsTable  = document.getElementById('region-stats-table');
const bcTable           = document.getElementById('bc-table');
const bcTableCaption    = document.getElementById('bc-table-caption');
const csvDownloadBtn    = document.getElementById('csv-download-btn');
const unmappedSection   = document.getElementById('unmapped-section');
const unmappedSummary   = document.getElementById('unmapped-summary');
const unmappedReadsEl   = document.getElementById('unmapped-reads');
const unmappedExportN   = document.getElementById('unmapped-export-n');
const fastqDownloadBtn  = document.getElementById('fastq-download-btn');
const errorMsg          = document.getElementById('error-msg');
const modal             = document.getElementById('sample-modal');
const modalTitle        = document.getElementById('modal-title');
const modalBody         = document.getElementById('modal-body');
const modalClose        = document.getElementById('modal-close');

const optAutodetect   = document.getElementById('opt-autodetect');
const optFlanking     = document.getElementById('opt-flanking');
const optBroad        = document.getElementById('opt-broad');
const optTrim         = document.getElementById('opt-trim');
const optTrimVal      = document.getElementById('opt-trim-val');
const optUnknownLens  = document.getElementById('opt-unknown-lens');
const optLimit        = document.getElementById('opt-limit');
const optLimitVal     = document.getElementById('opt-limit-val');

// ── Utilities ─────────────────────────────────────────────────────────────────

function rc(seq) {
  const comp = { A: 'T', T: 'A', C: 'G', G: 'C' };
  let out = '';
  for (let i = seq.length - 1; i >= 0; i--) out += comp[seq[i]] ?? seq[i];
  return out;
}

function escapeHtml(s) {
  return String(s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

function fmtNum(n) {
  return Number(n).toLocaleString();
}

function downloadBlob(content, filename, type) {
  const blob = new Blob([content], { type });
  const url  = URL.createObjectURL(blob);
  const a    = document.createElement('a');
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}

// ── Read display ──────────────────────────────────────────────────────────────

// Render a single read card as an HTML string, with barcodes highlighted.
// sample: {title, seq, qual, extractions: [{name, bcStr, bcStart, bcEnd, isRc}]}
function formatReadHTML(sample) {
  const { title, seq, extractions } = sample;
  const readId     = title.split(/\s/)[0];
  const displaySeq = escapeHtml(seq.slice(0, 1000));

  let html = `<div class="read-card">
    <span style="color:#999;">${escapeHtml(readId)}</span>
    <span style="color:#bbb;font-size:11px;margin-left:8px;">read length: ${seq.length}</span><br>
    <div class="seq-box">${displaySeq}</div>`;

  extractions.forEach((ext, i) => {
    const color    = BC_COLORS[i % BC_COLORS.length];
    const rcLabel  = ext.isRc ? ' <span style="color:#aaa;font-size:11px;">(RC)</span>' : '';
    const nameHtml = `<span style="color:${color};font-weight:bold;">${escapeHtml(ext.name)}</span>${rcLabel}: `;

    if (!ext.bcStr) {
      html += nameHtml + '<span style="color:#c0392b;">not found</span><br>';
    } else {
      const useSeq  = ext.isRc ? rc(seq) : seq;
      const ctxS    = Math.max(0, ext.bcStart - 12);
      const ctxE    = Math.min(useSeq.length, ext.bcEnd + 12);
      const pre     = escapeHtml(useSeq.slice(ctxS, ext.bcStart));
      const bc      = escapeHtml(ext.bcStr);
      const post    = escapeHtml(useSeq.slice(ext.bcEnd, ctxE));
      const lell    = ctxS > 0 ? '…' : '';
      const rell    = ctxE < useSeq.length ? '…' : '';
      const posLbl  = `<span style="color:#bbb;font-size:11px;margin-left:6px;">pos ${ext.bcStart}–${ext.bcEnd}</span>`;
      html += nameHtml
        + `${lell}<span style="color:#555;">${pre}</span>`
        + `<span style="background:${color};color:#fff;padding:1px 3px;border-radius:3px;font-weight:bold;">${bc}</span>`
        + `<span style="color:#555;">${post}</span>${rell}`
        + posLbl + '<br>';
    }
  });

  html += '</div>';
  return html;
}

// ── Parse options ─────────────────────────────────────────────────────────────

function getOpts() {
  const unknownLensVal = optUnknownLens.value;
  return {
    autodetect:  optAutodetect.checked,
    flankingLen: parseInt(optFlanking.value, 10) || 8,
    broadSearch: optBroad.checked,
    trimStart:   optTrim.checked ? (parseInt(optTrimVal.value, 10) || 0) : 0,
    unknownLens: unknownLensVal === 'all' ? 'All' : null,
    readCutoff:  optLimit.checked ? (parseInt(optLimitVal.value, 10) || null) : null,
    maxSamplesPerBc: 20,
    maxFailSamples:  1000,
    progressInterval: 100000,
  };
}

// ── Render results ────────────────────────────────────────────────────────────

function renderStatCards(stats, regions) {
  const { totalReads } = stats;
  const totalFails  = regions.reduce((s, r) => s + (stats[r.name]?.RegexFail ?? 0), 0);
  const nCombos     = currentResult.sortedCounts.filter(r =>
    regions.every(reg => r[reg.name] !== 'RegexFail')
  ).length;

  const cards = [
    { label: 'Total reads',        value: fmtNum(totalReads) },
    { label: 'Unique BC combos',   value: fmtNum(nCombos) },
    { label: 'Barcode regions',    value: regions.length },
    { label: 'Unmapped (stored)',  value: fmtNum(currentResult.failSamples.length) },
  ];

  statCards.innerHTML = cards.map(c =>
    `<div class="stat-card"><div class="label">${c.label}</div><div class="value">${c.value}</div></div>`
  ).join('');
}

function renderRegionStats(stats, regions) {
  if (!regions.length) { regionStatsSection.hidden = true; return; }
  const headerRow = `<tr><th>Region</th><th>Unique BCs</th><th>Regex fails</th></tr>`;
  const rows = regions.map(r => {
    const s = stats[r.name] ?? {};
    return `<tr><td>${escapeHtml(r.name)}</td><td>${fmtNum(s.nUnique ?? 0)}</td><td>${fmtNum(s.RegexFail ?? 0)}</td></tr>`;
  }).join('');
  regionStatsTable.innerHTML = `<thead>${headerRow}</thead><tbody>${rows}</tbody>`;
  regionStatsSection.hidden = false;
}

function renderBarcodeTable(sortedCounts, regions) {
  const top      = sortedCounts.slice(0, 1000);
  const colNames = regions.map(r => r.name).concat(['Count']);

  const headerCells = colNames.map(c => `<th>${escapeHtml(c)}</th>`).join('');
  const bodyRows    = top.map((row, i) => {
    const cells = colNames.map(c => `<td>${escapeHtml(String(row[c] ?? ''))}</td>`).join('');
    return `<tr data-row-index="${i}">${cells}</tr>`;
  }).join('');

  bcTable.innerHTML = `<thead><tr>${headerCells}</tr></thead><tbody>${bodyRows}</tbody>`;

  bcTableCaption.textContent =
    `Showing top ${fmtNum(top.length)} of ${fmtNum(sortedCounts.length)} unique combinations.`;

  // Row click → open modal
  bcTable.querySelectorAll('tbody tr').forEach(tr => {
    tr.addEventListener('click', () => {
      const idx  = parseInt(tr.dataset.rowIndex, 10);
      const row  = top[idx];
      const key  = regions.map(r => row[r.name]).join('\t');
      const samples = currentResult.readSamples.get(key) ?? [];
      const label   = regions.map(r => `${r.name}: ${row[r.name]}`).join('  |  ');
      openModal(label, samples);
    });
  });
}

function renderUnmapped(failSamples) {
  if (!failSamples.length) { unmappedSection.hidden = true; return; }

  unmappedSummary.textContent = `Unmapped reads (${fmtNum(failSamples.length)} stored)`;
  unmappedReadsEl.innerHTML   = failSamples.slice(0, 50).map(s => formatReadHTML(s)).join('');
  unmappedExportN.max         = failSamples.length;
  unmappedExportN.value       = Math.min(100, failSamples.length);
  unmappedSection.hidden      = false;

  fastqDownloadBtn.onclick = () => {
    const n    = Math.min(parseInt(unmappedExportN.value, 10) || 100, failSamples.length);
    const text = failSamples.slice(0, n)
      .map(s => `${s.title}\n${s.seq}\n+\n${s.qual}`).join('\n');
    downloadBlob(text, 'unmapped_reads.fastq', 'text/plain');
  };
}

function renderResults({ sortedCounts, readSamples, failSamples, stats, regions }) {
  currentResult = { sortedCounts, readSamples, failSamples, stats, regions };

  renderStatCards(stats, regions);
  renderRegionStats(stats, regions);
  renderBarcodeTable(sortedCounts, regions);
  renderUnmapped(failSamples);

  // CSV download
  csvDownloadBtn.onclick = () => {
    const cols   = regions.map(r => r.name).concat(['Count']);
    const header = cols.join(',');
    const body   = sortedCounts.map(row =>
      cols.map(c => JSON.stringify(String(row[c] ?? ''))).join(',')
    ).join('\n');
    downloadBlob(header + '\n' + body, 'barcode_counts.csv', 'text/csv');
  };

  errorMsg.hidden = true;
  resultsEl.hidden = false;
}

// ── Modal ─────────────────────────────────────────────────────────────────────

function openModal(label, samples) {
  modalTitle.textContent = label + (samples.length ? ` — ${samples.length} stored read(s)` : ' — no stored reads');
  modalBody.innerHTML    = samples.length
    ? samples.map(s => formatReadHTML(s)).join('')
    : '<p style="color:#888;font-size:13px;">No stored reads for this barcode.</p>';
  modal.showModal();
}

modalClose.addEventListener('click', () => modal.close());
modal.addEventListener('click', e => { if (e.target === modal) modal.close(); });

// ── Construct source toggle ───────────────────────────────────────────────────

document.querySelectorAll('input[name="construct-src"]').forEach(radio => {
  radio.addEventListener('change', () => {
    const isFasta = radio.value === 'fasta';
    constructTextWrap.hidden  = isFasta;
    constructFastaWrap.hidden = !isFasta;
    if (!isFasta) constructSeq = constructInput.value.trim() || BARSEQ_DEFAULT;
  });
});

constructInput.value = BARSEQ_DEFAULT;
constructInput.addEventListener('input', () => {
  constructSeq = constructInput.value.trim() || BARSEQ_DEFAULT;
});

fastaUpload.addEventListener('change', () => {
  const file = fastaUpload.files[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = e => {
    const content = e.target.result;
    constructSeq  = content.split('\n')
      .filter(l => !l.startsWith('>') && l.trim())
      .join('');
    fastaHint.textContent = `Loaded: ${file.name} (${constructSeq.length} bp)`;
  };
  reader.readAsText(file);
});

// ── Option toggles ────────────────────────────────────────────────────────────

optTrim.addEventListener('change', () => { optTrimVal.disabled = !optTrim.checked; });
optLimit.addEventListener('change', () => { optLimitVal.disabled = !optLimit.checked; });

// ── FASTQ upload → enable parse button ───────────────────────────────────────

fastqUpload.addEventListener('change', () => {
  const ready       = !!fastqUpload.files[0];
  parseBtn.disabled = !ready;
  parseHint.textContent = ready ? '' : 'Upload a FASTQ file to enable parsing.';
});

// ── Parse ─────────────────────────────────────────────────────────────────────

parseBtn.addEventListener('click', () => {
  const file = fastqUpload.files[0];
  if (!file) return;

  // Terminate any previous worker
  if (workerInstance) workerInstance.terminate();

  // Reset UI
  errorMsg.hidden     = true;
  resultsEl.hidden    = true;
  progressWrap.hidden = false;
  progressFill.style.width = '0%';
  progressLabel.textContent = 'Starting…';
  parseBtn.disabled = true;

  const isGzip = file.name.endsWith('.gz');
  const opts   = getOpts();

  file.arrayBuffer().then(bytes => {
    workerInstance = new Worker('js/worker.js');

    workerInstance.onmessage = ({ data }) => {
      if (data.type === 'progress') {
        progressLabel.textContent = `Parsed ${fmtNum(data.count)} reads…`;
        // Indeterminate progress — just animate
        progressFill.style.width = '60%';
      } else if (data.type === 'done') {
        progressWrap.hidden = true;
        parseBtn.disabled   = false;
        renderResults(data);
        workerInstance = null;
      } else if (data.type === 'error') {
        progressWrap.hidden     = true;
        parseBtn.disabled       = false;
        errorMsg.textContent    = 'Parsing failed: ' + data.message;
        errorMsg.hidden         = false;
        resultsEl.hidden        = false;
        workerInstance          = null;
      }
    };

    workerInstance.onerror = err => {
      progressWrap.hidden  = true;
      parseBtn.disabled    = false;
      errorMsg.textContent = 'Worker error: ' + err.message;
      errorMsg.hidden      = false;
      resultsEl.hidden     = false;
      workerInstance       = null;
    };

    workerInstance.postMessage({ bytes, isGzip, construct: constructSeq, opts });
  });
});
