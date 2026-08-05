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
  function highlightTriad(model, triad, isReference, entry) {
    const nums = [triad.ser, triad.asp, triad.his].filter(n => Number.isFinite(n));
    if (!nums.length) return;
    const radius = isReference ? 0.15 : 0.25;
    // add=true so the cartoon underneath survives.
    model.setStyle({ resi: nums }, { stick: { radius, colorscheme: 'default' } }, true);

    // A translucent surface over ALL THREE residues, not a sphere on the nucleophile.
    // The sphere marked one atom and said nothing about the shape of the pocket the three
    // residues form together, which is the thing being compared between structures.
    //
    // A surface belongs to the VIEWER, not to the model, so unlike setStyle it is not
    // undone by restyling and has to be tracked and removed explicitly. Its id is kept on
    // the entry so toggling the triad off, hiding a structure or removing it entirely all
    // clean up rather than leaving an orphaned surface floating in the scene.
    addTriadSurface(model, nums, entry);
  }

  /* Fluorescent yellow, chosen against the viewer's near-black background: the palette
     used for the structures themselves is all blues and warm mid-tones, so the triad needs
     a hue that belongs to no structure and reads at low opacity. */
  const TRIAD_SURFACE_COLOUR = '#eaff00';
  const TRIAD_SURFACE_OPACITY = 0.55;
  /* Mutated positions get their own surface in fluorescent pink. Two surfaces, two
     questions: yellow is "where the chemistry happens", pink is "what was changed", and
     the interesting thing about an engineered variant is how far the second sits from the
     first. Distinct enough from the yellow to read where a mutation is adjacent to the
     active site, which is exactly the case worth seeing. */
  const MUT_SURFACE_COLOUR = '#ff2ec4';
  const MUT_SURFACE_OPACITY = 0.5;

  function addTriadSurface(model, nums, entry) {
    if (!viewer || !entry) return;
    removeTriadSurface(entry);
    try {
      const sel = { resi: nums, model: model };
      const out = viewer.addSurface(
        $3Dmol.SurfaceType.VDW,
        { opacity: TRIAD_SURFACE_OPACITY, color: TRIAD_SURFACE_COLOUR },
        sel, sel);
      // 3Dmol changed addSurface from returning a surface id to returning a Promise that
      // resolves to one. Handle both so this does not silently stop cleaning up when the
      // CDN serves a newer build.
      if (out && typeof out.then === 'function') {
        out.then(res => {
          entry.surf = (res && res.surfid !== undefined) ? res.surfid : res;
          // MUST repaint. A surface is built asynchronously off the main thread and
          // nothing in 3Dmol repaints when it finishes, so a surface that completes after
          // the last render() is simply never drawn. This is why the triad appeared and
          // the mutations did not: the triad belonged to the reference, added first, and
          // a later add() happened to repaint it, while the surfaces created by the LAST
          // add() had no render after them at all.
          viewer.render();
        });
      } else {
        entry.surf = (out && out.surfid !== undefined) ? out.surfid : out;
        viewer.render();
      }
    } catch (err) {
      // A surface is an enhancement; failing to build one must not take the viewer down.
      entry.surf = null;
    }
  }

  function removeTriadSurface(entry) {
    if (!entry || entry.surf === null || entry.surf === undefined) return;
    try { viewer.removeSurface(entry.surf); } catch (err) { /* already gone */ }
    entry.surf = null;
  }

  /* Mutated residues, same mechanism as the triad surface and tracked separately so the
     two can be toggled independently. */
  function addMutationSurface(model, nums, entry) {
    if (!viewer || !entry || !nums || !nums.length) return;
    removeMutationSurface(entry);
    try {
      const sel = { resi: nums, model: model };
      const out = viewer.addSurface(
        $3Dmol.SurfaceType.VDW,
        { opacity: MUT_SURFACE_OPACITY, color: MUT_SURFACE_COLOUR },
        sel, sel);
      if (out && typeof out.then === 'function') {
        out.then(res => {
          entry.mutSurf = (res && res.surfid !== undefined) ? res.surfid : res;
          viewer.render();          // see the note in addTriadSurface
        });
      } else {
        entry.mutSurf = (out && out.surfid !== undefined) ? out.surfid : out;
        viewer.render();
      }
    } catch (err) {
      entry.mutSurf = null;
    }
  }

  function removeMutationSurface(entry) {
    if (!entry || entry.mutSurf === null || entry.mutSurf === undefined) return;
    try { viewer.removeSurface(entry.mutSurf); } catch (err) { /* already gone */ }
    entry.mutSurf = null;
  }

  /* Residue selection, shared with whatever else is on the page.
   *
   * The viewer and a sequence panel are two views of one molecule, so a click in either
   * should move the selection in both. This owns the viewer half: it makes the model's
   * atoms clickable, draws a marker on the picked residue, and reports the residue number
   * to a callback the page supplies.
   */
  let onPick = null;
  let pickedSurf = null;
  let pickedResi = null;
  const PICK_COLOUR = '#ffffff';

  function setResidueClickHandler(fn) { onPick = fn; }

  function makeClickable(model, id) {
    try {
      model.setClickable({}, true, (atom) => {
        if (!atom || !Number.isFinite(atom.resi)) return;
        // The id travels with the click: with several structures overlaid, "residue 160"
        // is ambiguous until you know whose.
        selectResidue(atom.resi, true, id);
      });
    } catch (err) { /* older builds without setClickable: clicking simply does nothing */ }
  }

  /** Mark a residue in the viewer. `notify` false when the call CAME from the page, so a
   *  sequence click does not bounce back and re-enter its own handler. */
  function selectResidue(resi, notify, id) {
    if (pickedSurf !== null && pickedSurf !== undefined) {
      try { viewer.removeSurface(pickedSurf); } catch (err) { /* gone */ }
      pickedSurf = null;
    }
    pickedResi = resi;
    if (Number.isFinite(resi)) {
      const target = (id && loaded.get(id)) ||
                     Array.from(loaded.values()).find(e => !e.isReference) ||
                     loaded.values().next().value;
      if (target) {
        const sel = { resi: [resi], model: target.model };
        try {
          const out = viewer.addSurface($3Dmol.SurfaceType.VDW,
            { opacity: 0.85, color: PICK_COLOUR }, sel, sel);
          if (out && typeof out.then === 'function') {
            out.then(res => {
              pickedSurf = (res && res.surfid !== undefined) ? res.surfid : res;
              viewer.render();
            });
          } else { pickedSurf = out; viewer.render(); }
        } catch (err) { /* selection marker is an enhancement */ }
        viewer.zoomTo(sel);
        viewer.zoom(0.6);
      }
    }
    viewer.render();
    if (notify && typeof onPick === 'function') onPick(resi, id);
  }

  function setMutationsVisible(on) {
    loaded.forEach(entry => {
      removeMutationSurface(entry);
      if (on && entry.mutations && entry.mutations.length) {
        addMutationSurface(entry.model, entry.mutations, entry);
      }
    });
    viewer.render();
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

    // The entry exists before the triad is drawn, because the surface id is stored on it.
    const entry = { model, colour, triad: opts.triad || null,
                    mutations: opts.mutations || null,
                    isReference: !!opts.isReference, surf: null, mutSurf: null };
    loaded.set(id, entry);
    if (opts.triad) highlightTriad(model, opts.triad, opts.isReference, entry);
    makeClickable(model, id);
    if (entry.mutations && entry.mutations.length) {
      // Sticks as well as the surface: the surface says where, the sticks say which side
      // chain, and at 0.5 opacity the sticks stay legible underneath.
      model.setStyle({ resi: entry.mutations },
                     { stick: { radius: 0.2, colorscheme: 'default' } }, true);
      addMutationSurface(model, entry.mutations, entry);
    }
    focusActiveSite();
    viewer.render();
    return id;
  }

  function setTriadVisible(on) {
    loaded.forEach(entry => {
      entry.model.setStyle({}, styleFor(entry.colour, entry.isReference));
      removeTriadSurface(entry);
      if (on && entry.triad) highlightTriad(entry.model, entry.triad, entry.isReference, entry);
    });
    viewer.render();
  }

  function remove(id) {
    const entry = loaded.get(id);
    if (!entry) return;
    removeTriadSurface(entry);
    removeMutationSurface(entry);
    viewer.removeModel(entry.model);
    loaded.delete(id);
    viewer.render();
  }

  function setVisible(id, visible) {
    const entry = loaded.get(id);
    if (!entry) return;
    entry.model.setStyle({}, visible ? styleFor(entry.colour, entry.isReference) : {});
    removeTriadSurface(entry);
    if (visible && entry.triad) highlightTriad(entry.model, entry.triad, entry.isReference, entry);
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

  /* ---- extra controls ---- */

  let spinning = false;
  function toggleSpin() {
    if (!viewer) return false;
    spinning = !spinning;
    // 3Dmol drives its own animation loop, so this does not need a rAF here.
    viewer.spin(spinning ? 'y' : false);
    return spinning;
  }

  /** Cartoon or sticks for the backbone, WITHOUT disturbing the highlights: the triad and
   *  mutation surfaces are separate objects, but their stick overlays are styles on the
   *  same model and would be wiped by a bare setStyle, so they are reapplied. */
  function setRepresentation(kind) {
    if (!viewer) return;
    loaded.forEach(entry => {
      const base = (kind === 'stick')
        ? { stick: { radius: entry.isReference ? 0.08 : 0.12, color: entry.colour } }
        : styleFor(entry.colour, entry.isReference);
      entry.model.setStyle({}, base);
      if (entry.triad) {
        const nums = [entry.triad.ser, entry.triad.asp, entry.triad.his].filter(Number.isFinite);
        if (nums.length) {
          entry.model.setStyle({ resi: nums },
            { stick: { radius: entry.isReference ? 0.15 : 0.25, colorscheme: 'default' } }, true);
        }
      }
      if (entry.mutations && entry.mutations.length) {
        entry.model.setStyle({ resi: entry.mutations },
          { stick: { radius: 0.2, colorscheme: 'default' } }, true);
      }
    });
    viewer.render();
  }

  /** Save the canvas as a PNG. Uses 3Dmol's own pngURI so the image is the rendered
   *  scene rather than a re-read of the framebuffer, which comes back blank on some
   *  drivers once the context has been presented. */
  function savePNG(filename) {
    if (!viewer) return;
    try {
      const uri = viewer.pngURI();
      const a = document.createElement('a');
      a.href = uri;
      a.download = filename || 'structure.png';
      document.body.appendChild(a); a.click(); a.remove();
    } catch (err) {
      showError('Could not save the image: ' + err.message);
    }
  }

  function reset() { if (viewer) { focusActiveSite(); viewer.render(); } }

  function zoomAll() { if (viewer) { viewer.zoomTo(); viewer.render(); } }

  function clear() { Array.from(loaded.keys()).forEach(remove); }

  return { init, add, remove, setVisible, setTriadVisible, setMutationsVisible,
           setResidueClickHandler, selectResidue,
           toggleSpin, setRepresentation, savePNG,
           reset, zoomAll, clear, PALETTE, REFERENCE_COLOUR };
})();
