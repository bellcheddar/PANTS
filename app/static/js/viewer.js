/* Superposed structure viewer, on 3Dmol.js.
 *
 * Two design points.
 *
 * 1. Every PDB served from /static/structures/ has ALREADY been superposed onto IsPETase
 *    (6EQE) at write time by pipeline/structure/fold.py. Overlaying candidates is just
 *    loading several files: no alignment happens in the browser. If structures ever look
 *    misaligned, the bug is in the pipeline, not here.
 *
 * 2. 3Dmol rather than Mol*. The Mol* UMD viewer bundle exports only `Viewer` (no Color,
 *    no MolScriptBuilder, no selection API), so per-structure colouring and residue
 *    selection are simply not reachable from it. 3Dmol exposes both directly, and is what
 *    the sibling AlphaFraud and BoltzMaker viewers already use.
 *
 * The cartoon depends on the HELIX/SHEET records the pipeline computes with biotite. This
 * 3Dmol build does not set atom.ss from them itself, so they are parsed here and applied.
 * Without that step every structure renders as a featureless tube.
 */

const PANTSViewer = (() => {
  'use strict';

  const PALETTE = ['#1e73be', '#ff6900', '#00d084', '#9b51e0', '#fcb900', '#e0245e'];
  const REFERENCE_COLOUR = '#9aa5b1';

  let viewer = null;
  const loaded = new Map();   // id -> { model, colour, triad, isReference }

  function showError(msg) {
    const el = document.getElementById('viewer-error');
    if (el) { el.textContent = msg; el.style.display = 'block'; }
    console.error('[PANTSViewer]', msg);
  }

  function init(elementId) {
    if (viewer) return viewer;
    const el = document.getElementById(elementId);
    if (!el) { showError('viewer container not found'); return null; }
    if (typeof $3Dmol === 'undefined') { showError('3Dmol failed to load'); return null; }
    viewer = $3Dmol.createViewer(el, { backgroundColor: '#101418' });
    return viewer;
  }

  /* HELIX/SHEET -> per-residue code. Fixed PDB columns rather than whitespace splitting:
     residue numbers can run into the chain id at wide numbering. */
  function ssMapFromPdb(pdb) {
    const map = {};
    pdb.split('\n').forEach(line => {
      if (line.startsWith('HELIX')) {
        const s = parseInt(line.substring(21, 25), 10), e = parseInt(line.substring(33, 37), 10);
        if (!isNaN(s) && !isNaN(e)) for (let i = s; i <= e; i++) map[i] = 'h';
      } else if (line.startsWith('SHEET')) {
        const s = parseInt(line.substring(22, 26), 10), e = parseInt(line.substring(33, 37), 10);
        if (!isNaN(s) && !isNaN(e)) for (let i = s; i <= e; i++) map[i] = 's';
      }
    });
    return map;
  }

  function applySS(model, pdb) {
    const map = ssMapFromPdb(pdb);
    model.selectedAtoms({}).forEach(a => { a.ss = map[a.resi] || 'c'; });
  }

  function styleFor(colour, isReference) {
    return isReference
      ? { cartoon: { color: colour, opacity: 0.5 } }
      : { cartoon: { color: colour } };
  }

  /**
   * Draw the catalytic triad as sticks, with a translucent sphere marking the nucleophile.
   *
   * Residue numbers come from the geometry stage, which found the triad by GEOMETRY (a
   * spatially connected Ser/His/Asp) rather than by sequence position, so what is drawn is
   * what was actually measured, not a guess from an alignment.
   */
  function highlightTriad(model, triad, isReference) {
    const nums = [triad.ser, triad.asp, triad.his].filter(n => Number.isFinite(n));
    if (!nums.length) return;
    const radius = isReference ? 0.15 : 0.25;
    // add=true so the cartoon underneath survives.
    model.setStyle({ resi: nums }, { stick: { radius, colorscheme: 'default' } }, true);
    if (Number.isFinite(triad.ser)) {
      model.setStyle({ resi: [triad.ser] }, {
        stick: { radius, colorscheme: 'default' },
        sphere: { radius: 0.8, color: '#ff2d55', opacity: 0.4 },
      }, true);
    }
  }

  /**
   * Load one structure.
   * @param {string} id     candidate id, or 'reference'
   * @param {string} url    PDB URL, already superposed
   * @param {object} opts   { colour, isReference, triad: {ser, asp, his} }
   */
  async function add(id, url, opts = {}) {
    if (loaded.has(id)) return id;
    if (!init('viewer')) return null;

    const colour = opts.colour || (opts.isReference
      ? REFERENCE_COLOUR : PALETTE[loaded.size % PALETTE.length]);

    let text;
    try {
      const resp = await fetch(url);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      text = await resp.text();
    } catch (err) {
      showError(`Could not fetch ${id}: ${err.message}`);
      return null;
    }

    const model = viewer.addModel(text, 'pdb');
    applySS(model, text);
    model.setStyle({}, styleFor(colour, opts.isReference));
    if (opts.triad) highlightTriad(model, opts.triad, opts.isReference);

    loaded.set(id, { model, colour, triad: opts.triad || null,
                     isReference: !!opts.isReference });
    focusActiveSite();
    viewer.render();
    return id;
  }

  function setTriadVisible(on) {
    loaded.forEach(entry => {
      entry.model.setStyle({}, styleFor(entry.colour, entry.isReference));
      if (on && entry.triad) highlightTriad(entry.model, entry.triad, entry.isReference);
    });
    viewer.render();
  }

  function remove(id) {
    const entry = loaded.get(id);
    if (!entry) return;
    viewer.removeModel(entry.model);
    loaded.delete(id);
    viewer.render();
  }

  function setVisible(id, visible) {
    const entry = loaded.get(id);
    if (!entry) return;
    entry.model.setStyle({}, visible ? styleFor(entry.colour, entry.isReference) : {});
    if (visible && entry.triad) highlightTriad(entry.model, entry.triad, entry.isReference);
    viewer.render();
  }

  /* Frame the active site rather than the whole model.
   *
   * zoomTo() with no selection fits everything, and ESMFold predictions carry long
   * disordered termini that are far from the fold: fitting those shrinks the actual
   * protein to a knot in the middle of the viewport. Framing the catalytic triad puts the
   * thing being compared at the centre, which is the point of a superposed view. */
  function focusActiveSite() {
    if (!viewer) return;
    const nums = [];
    loaded.forEach(e => {
      if (e.triad) [e.triad.ser, e.triad.asp, e.triad.his]
        .filter(n => Number.isFinite(n)).forEach(n => nums.push(n));
    });
    if (nums.length) viewer.zoomTo({ resi: nums });
    else viewer.zoomTo();
    viewer.zoom(0.55);          // pull back so the surrounding fold stays in frame
  }

  function reset() { if (viewer) { focusActiveSite(); viewer.render(); } }

  function zoomAll() { if (viewer) { viewer.zoomTo(); viewer.render(); } }

  function clear() { Array.from(loaded.keys()).forEach(remove); }

  return { init, add, remove, setVisible, setTriadVisible, reset, zoomAll, clear,
           PALETTE, REFERENCE_COLOUR };
})();
