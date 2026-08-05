/* Lineage figures for the PANTS home page.
 *
 * Plain SVG, no charting library: the droplet shares 3.9 GB with five other applications
 * and this page already carries 3Dmol for the structure viewers. Every view reads the same
 * payload the table below is rendered from, so a figure cannot drift from its own table.
 *
 * Values that do not exist stay missing rather than being imputed. Nine of the twenty-three
 * enzymes have no reported optimum temperature; each temperature view gives them a visible
 * "not reported" lane instead of dropping them, because a gap that is quietly removed reads
 * as a gap that is not there.
 */
(function () {
  const D = window.PANTS_FIGS;
  if (!D) return;
  const host = document.getElementById('fig-host');
  const cap = document.getElementById('fig-cap');
  const N = D.nodes;
  const TRIAD_C = '#eaff00', MUT_C = '#ff2ec4', ACC = '#8fb4ff', DIM = '#6e78a6';
  const esc = s => String(s).replace(/[&<>"]/g, c =>
    ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
  const el = (t, a, kids) => {
    const n = document.createElementNS('http://www.w3.org/2000/svg', t);
    for (const k in a) if (a[k] !== null && a[k] !== undefined) n.setAttribute(k, a[k]);
    (kids || []).forEach(c => n.appendChild(c));
    return n;
  };
  const txt = (x, y, s, cls, extra) => {
    const n = el('text', Object.assign({x, y, class: cls || 'fx-lab'}, extra || {}));
    n.textContent = s; return n;
  };
  const tip = (node, s) => { const t = el('title'); t.textContent = s; node.appendChild(t); return node; };
  const link = n => `/enzyme/${encodeURIComponent(n.id)}`;
  const go = n => { window.location.href = link(n); };
  function svg(w, h) { return el('svg', {viewBox: `0 0 ${w} ${h}`, role: 'img'}); }
  function mount(s, caption) { host.replaceChildren(s); cap.innerHTML = caption; }

  /* ---- 1. mutation hotspot track ---------------------------------------- */
  function hotspot() {
    const roots = Object.keys(D.hotspots).filter(r => Object.keys(D.hotspots[r]).length);
    const rows = [];
    roots.forEach(r => N.filter(n => n.root === r && n.muts.length)
      .sort((a, b) => a.nMut - b.nMut).forEach(n => rows.push(n)));
    const W = 940, L = 168, R = 24, rowH = 19, top = 34;
    const H = top + rows.length * rowH + 46;
    const s = svg(W, H);
    let y = top, lastRoot = null;
    roots.forEach(root => {
      const mine = rows.filter(n => n.root === root);
      const maxPos = Math.max(...mine.flatMap(n => n.muts), ...(D.triad[root] || [300]));
      const lo = 1, hi = Math.ceil((maxPos + 20) / 50) * 50;
      const X = p => L + (p - lo) / (hi - lo) * (W - L - R);
      s.appendChild(txt(4, y - 8, `${root} lineage · residue 1–${hi}`, 'fx-ttl'));
      (D.triad[root] || []).forEach(t => {
        s.appendChild(el('line', {x1: X(t), y1: y - 4, x2: X(t), y2: y + mine.length * rowH - 6,
          stroke: ACC, 'stroke-dasharray': '2 3', opacity: .45}));
      });
      mine.forEach(n => {
        s.appendChild(el('line', {x1: L, y1: y + 7, x2: W - R, y2: y + 7,
          stroke: '#26304d', opacity: .6}));
        const a = txt(4, y + 10, n.label.slice(0, 26), 'fx-lab fx-hit');
        a.style.cursor = 'pointer'; a.onclick = () => go(n); s.appendChild(a);
        n.muts.forEach(p => {
          const c = (D.hotspots[root] || {})[String(p)] || 1;
          const big = c >= 2;
          const t = el('line', {x1: X(p), y1: y + 1, x2: X(p), y2: y + 13,
            stroke: big ? TRIAD_C : DIM, 'stroke-width': big ? 2.4 : 1.2, class: 'fx-hit'});
          tip(t, `${n.label} · residue ${p}` + (big ? ` · shared by ${c} variants in this lineage` : ''));
          s.appendChild(t);
        });
        y += rowH;
      });
      Object.entries(D.hotspots[root] || {}).filter(([, c]) => c >= 3)
        .forEach(([p, c]) => s.appendChild(txt(X(+p) - 7, y + 11, p, 'fx-lab hi')));
      y += 30;
    });
    const top3 = Object.entries(D.hotspots.IsPETase || {}).sort((a, b) => b[1] - a[1]).slice(0, 4);
    mount(s, `One row per engineered variant, ticks at every substituted residue, along its own
      lineage's numbering. <b style="color:${TRIAD_C}">Yellow</b> marks a position two or more
      variants in that lineage share; blue dashes are the catalytic triad. The recurrence is the
      point: in the IsPETase lineage ${top3.map(([p, c]) => `<b>${p}</b> is hit by ${c}`).join(', ')}
      variants, from separate papers pursuing different goals. Click a name to open that enzyme.`);
  }

  /* ---- 2. pedigree ------------------------------------------------------ */
  function tree() {
    const W = 940, rowH = 22, pad = 26;
    const roots = [...new Set(N.map(n => n.root))]
      .sort((a, b) => N.filter(x => x.root === b).length - N.filter(x => x.root === a).length);
    const order = [];
    roots.forEach(r => {
      const fam = N.filter(n => n.root === r);
      const root = fam.find(n => n.depth === 0);
      if (root) order.push(root);
      fam.filter(n => n.depth === 1).sort((a, b) => (a.topt ?? 999) - (b.topt ?? 999))
        .forEach(n => { order.push(n); fam.filter(c => c.parent === n.id).forEach(c => order.push(c)); });
    });
    const H = pad + order.length * rowH + 20;
    const s = svg(W, H);
    const xOf = d => 30 + d * 46, yOf = i => pad + i * rowH;
    const idx = {}; order.forEach((n, i) => idx[n.id] = i);
    order.forEach((n, i) => {
      if (n.parent && idx[n.parent] !== undefined) {
        const py = yOf(idx[n.parent]), y = yOf(i);
        s.appendChild(el('path', {d: `M${xOf(n.depth - 1) + 5} ${py} L${xOf(n.depth - 1) + 5} ${y} L${xOf(n.depth) - 4} ${y}`,
          fill: 'none', stroke: '#26304d'}));
      }
      const warm = n.topt == null ? DIM
        : n.topt <= 40 ? '#00d084' : n.topt <= 55 ? '#fcb900' : '#ff4d5e';
      const c = el('circle', {cx: xOf(n.depth), cy: yOf(i), r: n.depth ? 4 : 5.5,
        fill: n.depth ? warm : 'none', stroke: warm, 'stroke-width': 1.6, class: 'fx-hit'});
      c.onclick = () => go(n);
      tip(c, `${n.id}${n.topt != null ? ' · ' + n.topt + ' °C' : ' · no optimum reported'}`);
      s.appendChild(c);
      const lab = txt(xOf(n.depth) + 11, yOf(i) + 3.5, n.label, 'fx-lab fx-hit');
      lab.style.cursor = 'pointer'; lab.onclick = () => go(n); s.appendChild(lab);
      if (n.nMut) s.appendChild(txt(xOf(n.depth) + 200, yOf(i) + 3.5, `+${n.nMut}`, 'fx-lab'));
      if (n.topt != null) s.appendChild(txt(xOf(n.depth) + 240, yOf(i) + 3.5, n.topt + ' °C', 'fx-lab'));
    });
    mount(s, `An engineering pedigree, <em>not</em> a phylogeny: every edge is a decision someone
      made, so branch length carries no evolutionary time. Node colour is optimum temperature —
      <b style="color:#00d084">≤40 °C</b>, <b style="color:#fcb900">≤55</b>,
      <b style="color:#ff4d5e">above</b>, grey where none is reported. Hollow rings are wild types.
      Note how many carry no engineering at all.`);
  }

  /* ---- 3. distance from 37 C ------------------------------------------- */
  function swarm() {
    const W = 940, H = 260, L = 40, R = 150, T = 40, B = 60;
    const withT = N.filter(n => n.topt != null), without = N.filter(n => n.topt == null);
    const lo = 25, hi = 82, X = t => L + (t - lo) / (hi - lo) * (W - L - R);
    const s = svg(W, H);
    s.appendChild(el('rect', {x: X(32), y: T - 12, width: X(42) - X(32), height: H - T - B + 18,
      fill: TRIAD_C, opacity: .07}));
    s.appendChild(el('line', {x1: X(37), y1: T - 14, x2: X(37), y2: H - B + 6,
      stroke: TRIAD_C, 'stroke-width': 1.4}));
    s.appendChild(txt(X(37) - 30, T - 20, '37 °C — the target', 'fx-lab hi'));
    s.appendChild(el('line', {x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: '#26304d'}));
    [30, 40, 50, 60, 70, 80].forEach(t =>
      s.appendChild(txt(X(t) - 7, H - B + 16, t + ' °C', 'fx-lab')));
    const used = {};
    withT.sort((a, b) => a.topt - b.topt).forEach(n => {
      const k = Math.round(X(n.topt) / 13);
      used[k] = (used[k] || 0) + 1;
      const y = H - B - 10 - (used[k] - 1) * 13;
      const near = Math.abs(n.topt - 37) <= 5;
      const c = el('circle', {cx: X(n.topt), cy: y, r: 5.5, fill: near ? TRIAD_C : ACC,
        opacity: near ? 1 : .62, stroke: n.src === 'pdb' ? 'none' : '#0d1120',
        'stroke-width': 1.2, class: 'fx-hit'});
      c.onclick = () => go(n);
      tip(c, `${n.id} · ${n.topt} °C · ${n.toptText || ''}`);
      s.appendChild(c);
    });
    s.appendChild(el('line', {x1: W - R + 14, y1: T - 10, x2: W - R + 14, y2: H - B,
      stroke: '#26304d', 'stroke-dasharray': '3 3'}));
    s.appendChild(txt(W - R + 22, T - 2, `no optimum reported (${without.length})`, 'fx-lab'));
    without.forEach((n, i) => {
      const c = el('circle', {cx: W - R + 30 + (i % 5) * 15, cy: T + 16 + Math.floor(i / 5) * 15,
        r: 4.5, fill: 'none', stroke: DIM, class: 'fx-hit'});
      c.onclick = () => go(n); tip(c, `${n.id} · no optimum reported`); s.appendChild(c);
    });
    mount(s, `Every enzyme by how far its optimum sits from body temperature. The shaded band is
      ±5 °C of the therapeutic target, and the emptiness around it is the reason this project
      exists. The ${without.length} enzymes with no reported optimum are shown at the right rather
      than dropped, because a gap made by quietly removing rows is not a gap.
      <b>Optima are lossy</b>: several of these were measured on different substrates under
      different assays — hover for the qualifying text.`);
  }

  /* ---- 4. cleft vs temperature ----------------------------------------- */
  function scatter() {
    const W = 940, H = 320, L = 56, R = 30, T = 26, B = 48;
    const pts = N.filter(n => n.cleft != null && n.topt != null);
    const X = c => L + (c - 7) / (27 - 7) * (W - L - R);
    const Y = t => H - B - (t - 25) / (82 - 25) * (H - T - B);
    const s = svg(W, H);
    s.appendChild(el('line', {x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: '#26304d'}));
    s.appendChild(el('line', {x1: L, y1: H - B, x2: L, y2: T, stroke: '#26304d'}));
    s.appendChild(el('line', {x1: L, y1: Y(37), x2: W - R, y2: Y(37), stroke: TRIAD_C,
      'stroke-dasharray': '4 4', opacity: .8}));
    s.appendChild(txt(W - R - 54, Y(37) - 5, '37 °C', 'fx-lab hi'));
    [10, 15, 20, 25].forEach(c => s.appendChild(txt(X(c) - 6, H - B + 16, c + ' Å', 'fx-lab')));
    [30, 45, 60, 75].forEach(t => s.appendChild(txt(6, Y(t) + 3, t + ' °C', 'fx-lab')));
    pts.forEach(n => {
      const exp = n.src === 'pdb';
      const c = el('circle', {cx: X(n.cleft), cy: Y(n.topt), r: 6,
        fill: exp ? n.colour : 'none', stroke: n.colour, 'stroke-width': 1.8,
        opacity: .9, class: 'fx-hit'});
      c.onclick = () => go(n);
      tip(c, `${n.id} · cleft ${n.cleft} Å from ${n.nArom} aromatics · ${n.topt} °C · `
        + (exp ? 'crystal structure' : 'prediction'));
      s.appendChild(c);
      if (n.cleft < 12 || n.topt >= 70) s.appendChild(txt(X(n.cleft) + 9, Y(n.topt) + 3, n.label, 'fx-lab'));
    });
    mount(s, `Cleft width against optimum temperature. <b>Filled</b> markers are crystal structures,
      <b>hollow</b> are predictions — the distinction is load-bearing, because the loops gating the
      cleft are exactly what folding models get least reliably. TurboPETase sits alone at the left:
      real geometry, but measured over three aromatic residues where the rest have five, so it is
      not strictly comparable. Colour is lineage; click any point.`);
  }

  /* ---- 5. parallel coordinates ----------------------------------------- */
  function para() {
    const W = 940, H = 300, T = 40, B = 46, L = 70, R = 60;
    const dims = [['topt', 'Topt °C'], ['cleft', 'Cleft Å'], ['pid', '%ID to WT'],
                  ['nMut', 'Mutations'], ['len', 'Length aa']];
    const ext = {};
    dims.forEach(([k]) => {
      const v = N.map(n => n[k]).filter(x => x != null);
      ext[k] = [Math.min(...v), Math.max(...v)];
    });
    const ax = i => L + i * ((W - L - R) / (dims.length - 1));
    const Y = (k, v) => H - B - (v - ext[k][0]) / ((ext[k][1] - ext[k][0]) || 1) * (H - T - B);
    const s = svg(W, H);
    dims.forEach(([k, lab], i) => {
      s.appendChild(el('line', {x1: ax(i), y1: T, x2: ax(i), y2: H - B, stroke: '#26304d'}));
      s.appendChild(txt(ax(i) - 22, T - 12, lab, 'fx-ttl'));
      s.appendChild(txt(ax(i) - 14, T - 1, String(Math.round(ext[k][1])), 'fx-lab'));
      s.appendChild(txt(ax(i) - 14, H - B + 12, String(Math.round(ext[k][0])), 'fx-lab'));
    });
    N.forEach(n => {
      const pts = dims.map(([k], i) => n[k] == null ? null : [ax(i), Y(k, n[k])]);
      if (pts.filter(Boolean).length < 2) return;
      const d = pts.filter(Boolean).map((p, i) => (i ? 'L' : 'M') + p[0] + ' ' + p[1]).join(' ');
      const path = el('path', {d, fill: 'none', stroke: n.colour, 'stroke-width': 1.4,
        opacity: .5, class: 'fx-hit'});
      path.onmouseenter = () => { path.setAttribute('stroke-width', 3); path.setAttribute('opacity', 1); };
      path.onmouseleave = () => { path.setAttribute('stroke-width', 1.4); path.setAttribute('opacity', .5); };
      path.onclick = () => go(n);
      tip(path, `${n.id}` + dims.map(([k, l]) => ` · ${l} ${n[k] ?? '—'}`).join(''));
      s.appendChild(path);
    });
    mount(s, `Every enzyme as one line across all five stored properties, coloured by lineage.
      A line breaks where a value is missing rather than being interpolated across the gap.
      Hover to isolate a profile: it is the fastest way to see an enzyme that is unusual on
      several axes at once rather than on any single one.`);
  }

  /* ---- 6. evidence gap -------------------------------------------------- */
  function sun() {
    const W = 940, H = 300, cx = 300, cy = 150;
    const s = svg(W, H);
    const roots = [...new Set(N.map(n => n.root))]
      .sort((a, b) => N.filter(x => x.root === b).length - N.filter(x => x.root === a).length);
    const tot = N.length; let a0 = -Math.PI / 2;
    const arc = (r0, r1, a, b, fill, op) => {
      const p = (r, t) => [cx + r * Math.cos(t), cy + r * Math.sin(t)];
      const [x0, y0] = p(r0, a), [x1, y1] = p(r0, b), [x2, y2] = p(r1, b), [x3, y3] = p(r1, a);
      const big = b - a > Math.PI ? 1 : 0;
      return el('path', {d: `M${x0} ${y0} A${r0} ${r0} 0 ${big} 1 ${x1} ${y1} L${x2} ${y2}
        A${r1} ${r1} 0 ${big} 0 ${x3} ${y3} Z`, fill, opacity: op, class: 'fx-hit'});
    };
    roots.forEach(root => {
      const fam = N.filter(n => n.root === root);
      const a1 = a0 + fam.length / tot * Math.PI * 2;
      const g = arc(46, 78, a0 + .006, a1 - .006, fam[0].colour, .55);
      tip(g, `${root} lineage · ${fam.length} enzymes`); s.appendChild(g);
      let b0 = a0;
      fam.forEach(n => {
        const b1 = b0 + 1 / tot * Math.PI * 2;
        const has = n.nPdb > 0;
        const w = arc(80, 80 + (has ? 12 + n.nPdb * 5 : 6), b0 + .004, b1 - .004,
          has ? ACC : '#26304d', has ? .95 : 1);
        w.onclick = () => go(n);
        tip(w, `${n.id} · ${n.nPdb} experimental structure${n.nPdb === 1 ? '' : 's'}`
          + (has ? '' : ' · predicted only'));
        s.appendChild(w);
        b0 = b1;
      });
      a0 = a1;
    });
    const noPdb = N.filter(n => !n.nPdb);
    s.appendChild(txt(560, 40, 'Enzymes with no experimental structure', 'fx-ttl'));
    noPdb.forEach((n, i) => {
      const t = txt(566, 60 + i * 15, '· ' + n.label.slice(0, 34), 'fx-lab fx-hit');
      t.style.cursor = 'pointer'; t.onclick = () => go(n); s.appendChild(t);
    });
    mount(s, `Inner ring is lineage, outer spoke length is how many crystal structures an enzyme
      has; grey stubs have none. The asymmetry is the finding: the deepest structural record
      belongs to the industrial, high-temperature enzymes, while
      <b>${noPdb.length} of ${N.length}</b> — including every gut-metagenome enzyme and several
      therapeutically interesting variants — rest on predictions alone. The field has
      structurally characterised what industry needed.`);
  }

  /* ---- 7. timeline ------------------------------------------------------ */
  function time() {
    const W = 940, H = 300, L = 52, R = 40, T = 30, B = 46;
    const pts = N.filter(n => n.year && n.topt != null);
    const ys = pts.map(n => n.year);
    const y0 = Math.min(...ys), y1 = Math.max(...ys);
    const X = y => L + (y - y0) / ((y1 - y0) || 1) * (W - L - R);
    const Y = t => H - B - (t - 25) / (82 - 25) * (H - T - B);
    const s = svg(W, H);
    s.appendChild(el('line', {x1: L, y1: H - B, x2: W - R, y2: H - B, stroke: '#26304d'}));
    s.appendChild(el('line', {x1: L, y1: Y(37), x2: W - R, y2: Y(37), stroke: TRIAD_C,
      'stroke-dasharray': '4 4', opacity: .8}));
    s.appendChild(txt(L + 4, Y(37) - 5, '37 °C', 'fx-lab hi'));
    for (let y = y0; y <= y1; y++) s.appendChild(txt(X(y) - 11, H - B + 16, String(y), 'fx-lab'));
    [30, 45, 60, 75].forEach(t => s.appendChild(txt(6, Y(t) + 3, t + ' °C', 'fx-lab')));
    pts.forEach(n => {
      const p = N.find(x => x.id === n.parent);
      if (p && p.year && p.topt != null) {
        s.appendChild(el('line', {x1: X(p.year), y1: Y(p.topt), x2: X(n.year), y2: Y(n.topt),
          stroke: n.colour, opacity: .35}));
      }
    });
    pts.forEach(n => {
      const c = el('circle', {cx: X(n.year), cy: Y(n.topt), r: 5, fill: n.colour,
        opacity: .9, class: 'fx-hit'});
      c.onclick = () => go(n); tip(c, `${n.id} · ${n.year} · ${n.topt} °C`);
      s.appendChild(c);
      s.appendChild(txt(X(n.year) + 8, Y(n.topt) + 3, n.label.slice(0, 18), 'fx-lab'));
    });
    mount(s, `Publication year against optimum temperature, with a line where one enzyme was built
      on another. Read it as several groups working separately rather than one march toward heat:
      <b>year is the publication date, not the date a design decision was made</b>, and several of
      these studies were optimising expression yield or substrate loading, with temperature simply
      coming along. Enzymes without a year or an optimum are absent from this view only.`);
  }

  const VIEWS = {hotspot, tree, swarm, scatter, para, sun, time};
  const tabs = document.getElementById('fig-tabs');

  /* The chosen view lives in the URL fragment, so a particular figure can be linked to
     directly -- useful in a write-up, and it means the back button behaves. */
  function show(name, push) {
    const v = VIEWS[name] ? name : 'hotspot';
    tabs.querySelectorAll('button').forEach(x =>
      x.setAttribute('aria-selected', String(x.dataset.v === v)));
    VIEWS[v]();
    if (push && location.hash !== '#fig-' + v) history.replaceState(null, '', '#fig-' + v);
  }
  tabs.addEventListener('click', e => {
    const b = e.target.closest('button[data-v]'); if (!b) return;
    show(b.dataset.v, true);
  });
  window.addEventListener('hashchange', () => show((location.hash || '').replace('#fig-', ''), false));
  show((location.hash || '').replace('#fig-', ''), false);
})();
