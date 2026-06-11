importScripts('https://cdn.jsdelivr.net/npm/fflate@0.8.2/umd/index.js');
importScripts('parser.js');

self.onmessage = ({ data: { bytes, isGzip, construct, opts } }) => {
  try {
    const raw  = isGzip ? fflate.gunzipSync(new Uint8Array(bytes)) : new Uint8Array(bytes);
    const text = new TextDecoder().decode(raw);

    const regions = parseConstruct(construct, opts);
    if (regions.length === 0) {
      self.postMessage({ type: 'error', message: 'No barcode regions found in construct. Check the construct sequence and options.' });
      return;
    }

    const result = runParse(text, regions, opts, (count) => {
      self.postMessage({ type: 'progress', count });
    });

    self.postMessage({ type: 'done', ...result });
  } catch (err) {
    self.postMessage({ type: 'error', message: err.message });
  }
};
