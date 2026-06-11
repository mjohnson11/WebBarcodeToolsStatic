// ── DNA utilities ────────────────────────────────────────────────────────────

function rc(seq) {
  const comp = { A: 'T', T: 'A', C: 'G', G: 'C' };
  let out = '';
  for (let i = seq.length - 1; i >= 0; i--) out += comp[seq[i]] ?? seq[i];
  return out;
}

// ── Construct parsing ─────────────────────────────────────────────────────────

// Returns an array of region objects:
//   {name, flanking_seq_left, flanking_seq_right, len}
// where len is a number or 'unknown_len'.
function parseConstruct(construct, { autodetect = false, flankingLen = 8, unknownLens = null } = {}) {
  const regions = [];

  if (autodetect) {
    // Split on runs of N, W (A/T), or S (C/G)
    const splits = construct.split(/([NWS]+)/);
    let regionIdx = 0;
    for (let i = 0; i < splits.length; i++) {
      if (/^[NWS]+$/.test(splits[i])) {
        regions.push({
          name: `unknown_region_${++regionIdx}`,
          flanking_seq_left:  (splits[i - 1] ?? '').slice(-flankingLen),
          flanking_seq_right: (splits[i + 1] ?? '').slice(0, flankingLen),
          len: splits[i].length,
        });
      }
    }
  } else {
    // Named regions: (name:SEQUENCE)
    // Split with capturing groups yields: [before, name, seq, between, name, seq, ..., after]
    const splits = construct.split(/\((.*?):(.*?)\)/);
    const nRegions = Math.floor((splits.length - 1) / 3);
    for (let i = 0; i < nRegions; i++) {
      const before = splits[i * 3];
      const name   = splits[i * 3 + 1];
      const seq    = splits[i * 3 + 2];
      const after  = splits[(i + 1) * 3];
      regions.push({
        name,
        flanking_seq_left:  before.slice(-flankingLen),
        flanking_seq_right: after.slice(0, flankingLen),
        len: seq.length,
      });
    }
  }

  if (unknownLens === 'All') {
    regions.forEach(r => { r.len = 'unknown_len'; });
  } else if (Array.isArray(unknownLens)) {
    regions.forEach(r => { if (unknownLens.includes(r.name)) r.len = 'unknown_len'; });
  }

  return regions;
}

// ── Barcode extraction ────────────────────────────────────────────────────────

// Try to find exactly: leftFlank + bcLen chars + rightFlank in seq.
// Returns {bcStart, bcEnd, bcStr} or null.
function tryExtract(seq, leftFlank, rightFlank, bcLen) {
  let pos = seq.indexOf(leftFlank);
  while (pos !== -1) {
    const bcStart = pos + leftFlank.length;
    const bcEnd   = bcStart + bcLen;
    if (seq.slice(bcEnd, bcEnd + rightFlank.length) === rightFlank) {
      return { bcStart, bcEnd, bcStr: seq.slice(bcStart, bcEnd) };
    }
    pos = seq.indexOf(leftFlank, pos + 1);
  }
  return null;
}

// Extract barcode from seq (forward direction only).
// For known-length: cascades exact → ±1 → ±2 length with exact flanks.
// For unknown_len: greedy capture between exact left and right flanks.
function extractBarcode(seq, region) {
  const { flanking_seq_left: lf, flanking_seq_right: rf, len } = region;

  if (len === 'unknown_len') {
    let pos = seq.indexOf(lf);
    while (pos !== -1) {
      const bcStart = pos + lf.length;
      const rfPos   = seq.indexOf(rf, bcStart);
      if (rfPos !== -1) {
        return { bcStart, bcEnd: rfPos, bcStr: seq.slice(bcStart, rfPos) };
      }
      pos = seq.indexOf(lf, pos + 1);
    }
    return null;
  }

  for (const delta of [0, -1, 1, -2, 2]) {
    const tryLen = len + delta;
    if (tryLen < 1) continue;
    const result = tryExtract(seq, lf, rf, tryLen);
    if (result) return result;
  }
  return null;
}

// Try forward, then RC. Returns {bcStart, bcEnd, bcStr, isRc} or null.
function extractBarcodeWithRC(seq, region, broadSearch) {
  const fwd = extractBarcode(seq, region);
  if (fwd) return { ...fwd, isRc: false };
  if (broadSearch) {
    const rcSeq = rc(seq);
    const rev   = extractBarcode(rcSeq, region);
    if (rev) return { ...rev, isRc: true };
  }
  return null;
}

// ── FASTQ parsing ─────────────────────────────────────────────────────────────

function* parseFastq(text) {
  let start = 0;
  while (start < text.length) {
    let e1 = text.indexOf('\n', start);          if (e1 === -1) break;
    let e2 = text.indexOf('\n', e1 + 1);         if (e2 === -1) break;
    let e3 = text.indexOf('\n', e2 + 1);         if (e3 === -1) break;
    let e4 = text.indexOf('\n', e3 + 1);
    if (e4 === -1) e4 = text.length;

    const title = text.slice(start, e1).trimEnd();
    const seq   = text.slice(e1 + 1, e2).trimEnd();
    const plus  = text.slice(e2 + 1, e3).trimEnd();
    const qual  = text.slice(e3 + 1, e4).trimEnd();

    if (title.startsWith('@') && plus.startsWith('+') && seq.length > 0 && seq.length === qual.length) {
      yield { title, seq, qual };
    }
    start = e4 + 1;
  }
}

// ── Main parse loop ───────────────────────────────────────────────────────────

function runParse(text, regions, opts = {}, onProgress = null) {
  const {
    broadSearch      = true,
    trimStart        = 0,
    readCutoff       = null,
    maxSamplesPerBc  = 20,
    maxFailSamples   = 1000,
    progressInterval = 100000,
  } = opts;

  const counts         = new Map();  // bcKey (tab-joined) -> count
  const readSamples    = new Map();  // bcKey -> [{title, seq, qual, extractions}]
  const failSamples    = [];
  const regionUniqBcs  = regions.map(() => new Set());
  const regionFailCts  = regions.map(() => 0);

  let totalReads = 0;

  for (const { title, seq: rawSeq, qual: rawQual } of parseFastq(text)) {
    if (readCutoff !== null && totalReads >= readCutoff) break;

    const seq  = trimStart ? rawSeq.slice(trimStart) : rawSeq;
    const qual = trimStart ? rawQual.slice(trimStart) : rawQual;

    const extractions = regions.map((region, i) => {
      const hit = extractBarcodeWithRC(seq, region, broadSearch);
      if (hit) {
        regionUniqBcs[i].add(hit.bcStr);
        return { name: region.name, bcStr: hit.bcStr, bcStart: hit.bcStart, bcEnd: hit.bcEnd, isRc: hit.isRc };
      }
      regionFailCts[i]++;
      return { name: region.name, bcStr: null, bcStart: null, bcEnd: null, isRc: false };
    });

    // Use TAB as separator — can't appear in DNA sequences
    const bcKey = extractions.map(e => e.bcStr ?? 'RegexFail').join('\t');
    counts.set(bcKey, (counts.get(bcKey) ?? 0) + 1);

    const isFail = extractions.some(e => !e.bcStr);
    if (isFail) {
      if (failSamples.length < maxFailSamples) {
        failSamples.push({ title, seq, qual, extractions });
      }
    } else {
      if (!readSamples.has(bcKey)) readSamples.set(bcKey, []);
      const arr = readSamples.get(bcKey);
      if (arr.length < maxSamplesPerBc) {
        arr.push({ title, seq, qual, extractions });
      }
    }

    totalReads++;
    if (onProgress && totalReads % progressInterval === 0) onProgress(totalReads);
  }

  const stats = { totalReads };
  regions.forEach((r, i) => {
    stats[r.name] = { nUnique: regionUniqBcs[i].size, RegexFail: regionFailCts[i] };
  });

  const sortedCounts = [...counts.entries()]
    .sort((a, b) => b[1] - a[1])
    .map(([key, count]) => {
      const parts = key.split('\t');
      const row   = { Count: count };
      regions.forEach((r, i) => { row[r.name] = parts[i]; });
      return row;
    });

  return { sortedCounts, readSamples, failSamples, stats, regions };
}
