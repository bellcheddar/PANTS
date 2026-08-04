/* Superposed structure viewer, built on Mol*.
 *
 * The design point that makes this cheap: every mmCIF served from /static/structures/ has
 * ALREADY been superposed onto IsPETase (6EQE) at write time by pipeline/structure/fold.py.
 * So overlaying candidates is just loading several files into one viewer. No alignment
 * happens in the browser, which is what keeps an N-structure overlay responsive on a
 * phone and avoids shipping an alignment library to the client.
 *
 * If a structure ever appears misaligned, the bug is in the pipeline, not here.
 */

const PANTSViewer = (() => {
  'use strict';

  // Distinct, colour-blind-safe-ish palette. The reference is deliberately grey so
  // candidates read as the subject and IsPETase as the backdrop.
  const PALETTE = ['#1e73be', '#ff6900', '#00d084', '#9b51e0', '#fcb900', '#e0245e'];
  const REFERENCE_COLOUR = '#9aa5b1';

  let plugin = null;
  const loaded = new Map();   // id -> { ref, colour, triad }

  async function init(elementId) {
    if (plugin) return plugin;
    plugin = await molstar.Viewer.create(elementId, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: false,
      layoutShowLog: false,
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      pdbProvider: 'rcsb',
    });
    return plugin;
  }

  function colourFor(index, isReference) {
    return isReference ? REFERENCE_COLOUR : PALETTE[index % PALETTE.length];
  }

  /**
   * Load one structure.
   * @param {string} id        candidate id (or 'reference')
   * @param {string} url       mmCIF/PDB URL, already superposed
   * @param {object} opts      { colour, isReference, triad: {ser, asp, his} }
   */
  function showError(msg) {
    const el = document.getElementById('viewer-error');
    if (el) { el.textContent = msg; el.style.display = 'block'; }
    console.error('[PANTSViewer]', msg);
  }

  async function add(id, url, opts = {}) {
    if (loaded.has(id)) return;
    await init('viewer');
    try {

    const format = url.endsWith('.pdb') ? 'pdb' : 'mmcif';
    const colour = opts.colour || colourFor(loaded.size, opts.isReference);

    const data = await plugin.plugin.builders.data.download({ url, isBinary: false });
    const trajectory = await plugin.plugin.builders.structure.parseTrajectory(data, format);
    const model = await plugin.plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.plugin.builders.structure.createStructure(model);

    // applyPreset rather than addRepresentation: the preset builds the polymer component
    // first, which is what makes a cartoon possible. Adding a bare 'cartoon'
    // representation to the whole structure renders every atom as lines instead, because
    // nothing has classified the chains as polymer yet.
    const value = parseInt(colour.slice(1), 16);
    const cartoon = await plugin.plugin.builders.structure.representation.applyPreset(
      structure, 'polymer-cartoon',
      { theme: { globalName: 'uniform', globalColorParams: { value } } }
    );

    loaded.set(id, { structure, cartoon, colour, triad: opts.triad || null, reps: [] });

    if (opts.triad) await highlightTriad(id, opts.triad);
    return id;
    } catch (err) {
      showError(`Could not load ${id}: ${err && err.message ? err.message : err}`);
      throw err;
    }
  }

  /**
   * Draw the catalytic triad as ball-and-stick.
   *
   * Residue numbers come from the geometry table, which found them by GEOMETRY (a
   * spatially connected Ser/His/Asp) rather than by sequence position, so what is drawn
   * here is what was actually measured, not a guess from an alignment.
   */
  async function highlightTriad(id, triad) {
    const entry = loaded.get(id);
    if (!entry || !triad) return;
    const nums = [triad.ser, triad.asp, triad.his].filter(n => Number.isFinite(n));
    if (!nums.length) return;
    const MS = molstar.MolScriptBuilder;
    if (!MS) {
      console.warn('[PANTSViewer] MolScriptBuilder unavailable; triad not drawn');
      return;
    }

    const expr = MS.struct.generator.atomGroups({
      'residue-test': MS.core.set.has([
        MS.set(...nums),
        MS.struct.atomProperty.macromolecular.auth_seq_id(),
      ]),
    });

    const sel = await plugin.plugin.builders.structure.tryCreateComponentFromExpression(
      entry.structure, expr, `triad-${id}`, { label: `Catalytic triad (${id})` }
    );
    if (sel) {
      const rep = await plugin.plugin.builders.structure.representation.addRepresentation(sel, {
        type: 'ball-and-stick',
        color: 'element-symbol',
        typeParams: { sizeFactor: 0.35 },
      });
      entry.reps.push(rep);
    }
  }

  async function remove(id) {
    const entry = loaded.get(id);
    if (!entry) return;
    await plugin.plugin.build().delete(entry.structure).commit();
    loaded.delete(id);
  }

  async function setVisible(id, visible) {
    const entry = loaded.get(id);
    if (!entry) return;
    molstar.setSubtreeVisibility(plugin.plugin.state.data, entry.structure.ref, !visible);
  }

  function reset() { if (plugin) plugin.plugin.managers.camera.reset(); }

  async function clear() {
    for (const id of Array.from(loaded.keys())) await remove(id);
  }

  return { init, add, remove, setVisible, reset, clear, colourFor, PALETTE, REFERENCE_COLOUR };
})();
