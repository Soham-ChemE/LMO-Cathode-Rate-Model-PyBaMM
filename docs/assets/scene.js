/* Hero scene for the LMO rate-model page: the half-cell sandwich with a lithium-ion stream, the four measured discharge
   curves as 3D ribbons that draw as the reader scrolls (DFN model dashed through each), and the 2C loss split growing at
   the end of the track. Data: window.LMO from assets/data.js (built from the raw EC-Lab files and the fitted models). */
import { THREE, createStage, cameraPath } from './three-kit.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
const D = window.LMO, S = window.SITE, reduced = S.reduced, phone = S.phone;
const RATE_COL = { 'C/10': 0x22d3ee, 'C/5': 0x34d399, 'C/2': 0xf59e0b, '2C': 0xf87171 }, RATE_SPEED = { 'C/10': .12, 'C/5': .25, 'C/2': .5, '2C': 1.4 };
const textSprite = (txt, color, w = 2.4) => { const c = document.createElement('canvas'); c.width = 512; c.height = 96; const g = c.getContext('2d'); g.font = '600 40px -apple-system,Segoe UI,Roboto,Helvetica,Arial'; g.fillStyle = color; g.textAlign = 'center'; g.fillText(txt, 256, 62); const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), transparent: true, depthWrite: false })); s.scale.set(w, w * 96 / 512, 1); return s; };

(function () {
  const canvas = document.getElementById('heroCanvas'); if (!canvas) return;
  const st = createStage(canvas, { fov: 40, bloom: { strength: .35, radius: .45, threshold: .85 }, dof: { aperture: .0004, maxblur: .003 }, exposure: 1.0, envIntensity: .45 });
  if (st.failed) return; const { scene, camera } = st;
  const ctr = new OrbitControls(camera, canvas); ctr.enableDamping = true; ctr.enableZoom = false; ctr.enablePan = false; ctr.autoRotate = !reduced; ctr.autoRotateSpeed = .4; ctr.minPolarAngle = .85; ctr.maxPolarAngle = 1.5;
  const hemi = new THREE.HemisphereLight(0xffffff, 0x334466, .9); scene.add(hemi); const dl = new THREE.DirectionalLight(0xffffff, 1.5); dl.position.set(4, 8, 6); scene.add(dl);
  const rnd = (() => { let s = 7; return () => { s = (s * 16807) % 2147483647; return (s - 1) / 2147483646; }; })();
  // ---------- cell sandwich ----------
  const cell = new THREE.Group(); scene.add(cell);
  const H = 3.2, Dz = 2.2; const slab = (w, color, x, opts = {}) => { const m = new THREE.Mesh(new THREE.BoxGeometry(w, H, Dz), new THREE.MeshStandardMaterial({ color, roughness: .4, metalness: opts.metal || 0, transparent: !!opts.op, opacity: opts.op || 1 })); m.position.x = x; cell.add(m); return m; };
  slab(.5, 0xb8c0cc, -2.0, { metal: .7 }); slab(.25, 0xffffff, -1.55, { op: .28 }); slab(.12, 0x9aa5b4, .95, { metal: .85 });
  const lab = (t, x, y) => { const s = textSprite(t, '#8b98a9', 1.6); s.position.set(x, y, Dz / 2 + .2); cell.add(s); };
  lab('Li metal', -2, -H / 2 - .3); lab('separator', -1.55, H / 2 + .3); lab(`LiMn₂O₄ coating, ${D.cell ? D.cell.coating_um : 85} µm`, -.25, -H / 2 - .3); lab('Al', .95, H / 2 + .3);
  const NP = phone ? 160 : 260, pg = new THREE.SphereGeometry(.11, 12, 10), pm = new THREE.MeshStandardMaterial({ color: 0x7c5cff, roughness: .5 }); const parts = new THREE.InstancedMesh(pg, pm, NP); const PP = []; const m4 = new THREE.Matrix4();
  for (let i = 0; i < NP; i++) { const p = new THREE.Vector3(-1.35 + rnd() * 2.2, -H / 2 + .15 + rnd() * (H - .3), -Dz / 2 + .15 + rnd() * (Dz - .3)); PP.push(p); const sc = .75 + rnd() * .6; m4.makeScale(sc, sc, sc).setPosition(p); parts.setMatrixAt(i, m4); } cell.add(parts);
  const ely = new THREE.Mesh(new THREE.BoxGeometry(2.3, H, Dz), new THREE.MeshStandardMaterial({ color: 0x22d3ee, transparent: true, opacity: .08, depthWrite: false })); ely.position.x = -.25; cell.add(ely);
  const NI = phone ? 40 : 90, ig = new THREE.SphereGeometry(.05, 8, 6), im = new THREE.MeshStandardMaterial({ color: 0x22d3ee, emissive: 0x22d3ee, emissiveIntensity: 1.6 }); const ions = new THREE.InstancedMesh(ig, im, NI); cell.add(ions);
  const ION = []; for (let i = 0; i < NI; i++) ION.push({ t: rnd(), y: -H / 2 + .2 + rnd() * (H - .4), z: -Dz / 2 + .2 + rnd() * (Dz - .4), target: PP[Math.floor(rnd() * NP)] });
  let rate = 'C/2', dir = 1, dirT = 0; const rateBtns = [...document.querySelectorAll('.rates button')];
  rateBtns.forEach(b => b.addEventListener('click', () => { rate = b.dataset.rate; rateBtns.forEach(x => x.classList.toggle('on', x === b)); highlight(); st.dirty = true; }));
  function stepIons(dt) { const v = RATE_SPEED[rate] * dt; for (let i = 0; i < NI; i++) { const o = ION[i]; o.t += v * (.6 + .8 * ((i * 7) % 10) / 10); if (o.t > 1) { o.t = 0; o.target = PP[Math.floor(rnd() * NP)]; o.y = -H / 2 + .2 + rnd() * (H - .4); o.z = -Dz / 2 + .2 + rnd() * (Dz - .4); } const u = dir > 0 ? o.t : 1 - o.t; const x0 = -1.74, x1 = o.target.x; const x = x0 + (x1 - x0) * u; const y = o.y + (o.target.y - o.y) * Math.max(0, (u - .35) / .65); const z = o.z + (o.target.z - o.z) * Math.max(0, (u - .35) / .65); m4.makeTranslation(x, y, z); ions.setMatrixAt(i, m4); } ions.instanceMatrix.needsUpdate = true; }
  // ---------- voltage ribbons that draw with scroll ----------
  const rib = new THREE.Group(); scene.add(rib);
  const X = v => v / 1.5 * 5 - 2.5, Y = v => (v - 3.4) / 1.2 * 3, Z = i => 1.5 - i * 1.0; const ribbons = {};
  D.curves.forEach((c, i) => {
    const pts = c.measured; const pos = [], idx = []; for (let k = 0; k < pts.length; k++) { const x = X(pts[k][0]), y = Y(pts[k][1]), z = Z(i); pos.push(x, y, z - .14, x, y, z + .14); } for (let k = 0; k < pts.length - 1; k++) { const a = 2 * k; idx.push(a, a + 1, a + 2, a + 1, a + 3, a + 2); }
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3)); g.setIndex(idx); g.computeVertexNormals();
    const mesh = new THREE.Mesh(g, new THREE.MeshStandardMaterial({ color: RATE_COL[c.rate], side: THREE.DoubleSide, transparent: true, opacity: .92, roughness: .35, emissive: RATE_COL[c.rate], emissiveIntensity: .25 })); rib.add(mesh);
    const mp = c.model.map(p => new THREE.Vector3(X(p[0]), Y(p[1]) + .02, Z(i))); const line = new THREE.Line(new THREE.BufferGeometry().setFromPoints(mp), new THREE.LineDashedMaterial({ color: 0xffffff, dashSize: .12, gapSize: .08, transparent: true, opacity: .9 })); line.computeLineDistances(); rib.add(line);
    const s = textSprite(c.rate + (c.role === 'blind prediction' ? ' (blind)' : c.role === 'training' ? ' (fit)' : ''), '#' + RATE_COL[c.rate].toString(16).padStart(6, '0'), 1.5); s.position.set(X(pts[pts.length - 1][0]) + .7, Y(pts[pts.length - 1][1]), Z(i)); rib.add(s);
    ribbons[c.rate] = { mesh, line, label: s, nIdx: idx.length, nLine: mp.length };
  });
  const ax = new THREE.Group(); const lm = new THREE.LineBasicMaterial({ color: 0x8b98a9, transparent: true, opacity: .5 });
  const P = (a, b, c) => new THREE.Vector3(a, b, c); [[P(X(0), Y(3.5), Z(3) - .5), P(X(1.5), Y(3.5), Z(3) - .5)], [P(X(0), Y(3.5), Z(3) - .5), P(X(0), Y(4.5), Z(3) - .5)], [P(X(0), Y(3.5), Z(3) - .5), P(X(0), Y(3.5), Z(0) + .5)]].forEach(([a, b]) => ax.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints([a, b]), lm)));
  const axl = [textSprite('capacity 0 → 1.5 mAh', '#8b98a9', 2.2), textSprite('3.5 → 4.5 V', '#8b98a9', 1.4)]; axl[0].position.set(X(.75), Y(3.5) - .3, Z(3) - .5); axl[1].position.set(X(0) - .9, Y(4.0), Z(3) - .5); axl.forEach(s => ax.add(s)); rib.add(ax);
  function highlight() { for (const r in ribbons) { const on = r === rate; ribbons[r].mesh.material.opacity = on ? 1 : .5; ribbons[r].mesh.material.emissiveIntensity = on ? .6 : .15; ribbons[r].line.material.opacity = on ? 1 : .4; } } highlight();
  st.onTheme = dark => { hemi.intensity = .9 * st.lights; dl.intensity = 1.5 * st.lights; st.dirty = true; }; st.onTheme(st.dark);
  // ---------- 2C loss split: four bars that grow at the end of the track ----------
  const sh = D.stats.shares; const SPLIT = { o: Math.round(sh.contact_ohmic), r: Math.round(sh.reaction), d: Math.round(sh.particle_diffusion), e: Math.round(sh.electrolyte_concentration + sh.electrolyte_ohmic + sh.solid_ohmic) };
  const bars = new THREE.Group(); rib.add(bars); const BARC = { o: 0xf87171, r: 0xf59e0b, d: 0x34d399, e: 0x22d3ee }; const BARL = { o: 'ohmic + contact', r: 'charge transfer', d: 'particle diffusion', e: 'electrolyte' };
  const barMeshes = []; let bx = X(0) + .3;
  for (const k of ['o', 'r', 'd', 'e']) { const h = SPLIT[k] / 100 * 2.6; const m = new THREE.Mesh(new THREE.BoxGeometry(.55, h, .5), new THREE.MeshStandardMaterial({ color: BARC[k], emissive: BARC[k], emissiveIntensity: .3, roughness: .4 })); m.position.set(bx, Y(3.5) + h / 2, Z(3) - 1.6); m.scale.y = .001; const l = textSprite(`${SPLIT[k]} % ${BARL[k]}`, '#' + BARC[k].toString(16).padStart(6, '0'), 1.7); l.position.set(bx, Y(3.5) + h + .25, Z(3) - 1.6); l.visible = false; bars.add(m, l); barMeshes.push({ m, l, h }); bx += 1.25; }
  const barTitle = textSprite(`2C mean loss ${D.stats.two_c_mean_loss_mV} mV, split by the DFN`, '#e6edf3', 3.2); barTitle.position.set(X(0) + 2.2, Y(3.5) + 3.2, Z(3) - 1.6); barTitle.visible = false; bars.add(barTitle);
  // ---------- layout and camera path ----------
  const wide = () => st.wide !== false;
  function layout() { const w = wide(); cell.position.set(w ? -3.9 : -1.1, w ? -.2 : 1.4, 0); cell.scale.setScalar(w ? 1 : .55); rib.position.set(w ? 2.3 : 0, w ? -.6 : -2.2, 0); rib.scale.setScalar(w ? 1 : .6); }
  const pathWide = cameraPath([[-0.8, 2.6, 12.8], [2.5, 3.0, 10.4], [4.0, 1.6, 8.6], [3.6, 2.4, 7.4]], [[-1.4, -.5, 0], [.5, -.8, 0], [2.3, -.6, 0], [2.6, -.2, -.6]]);
  const pathNarrow = cameraPath([[0, 4.5, 14], [1.5, 3.5, 12], [0, 1.0, 9.5], [0, 1.5, 8.5]], [[0, .4, 0], [0, -.5, 0], [0, -1.6, 0], [0, -1.2, -.5]]);
  const pos = new THREE.Vector3(), tgt = new THREE.Vector3(); let dragging = false, lastT = -1, last = performance.now();
  canvas.addEventListener('pointerdown', () => dragging = true); addEventListener('pointerup', () => setTimeout(() => dragging = false, 2500)); ctr.addEventListener('change', () => st.dirty = true);
  layout(); (wide() ? pathWide : pathNarrow).at(reduced ? 1 : 0, pos, tgt); camera.position.copy(pos); ctr.target.copy(tgt); ctr.update();
  st.loop((t, dtMs) => {
    const dt = Math.min(.05, dtMs / 1000); const T = S.hero.t; let changed = false; layout();
    (wide() ? pathWide : pathNarrow).at(reduced ? 1 : T, pos, tgt);
    if (!dragging && camera.position.distanceToSquared(pos) > 1e-6) { camera.position.lerp(pos, reduced ? 1 : .08); changed = true; }
    ctr.target.lerp(tgt, .08); ctr.autoRotate = !reduced && T < .08 && !dragging; if (ctr.autoRotate || dragging) changed = true; ctr.update();
    if (!reduced) { dirT += dt; if (dirT > 7) { dirT = 0; dir *= -1; } stepIons(dt); changed = true; } else if (!ions.userData.once) { ions.userData.once = 1; stepIons(0); }
    // ribbons draw between 10 % and 65 % of the track; bars grow between 70 % and 95 %
    const draw = reduced ? 1 : Math.max(0, Math.min(1, (T - .1) / .55)); let i = 0;
    for (const r in ribbons) { const rb = ribbons[r]; const lead = Math.max(0, Math.min(1, draw * 1.25 - i * .08)); rb.mesh.geometry.setDrawRange(0, Math.floor(rb.nIdx * lead)); rb.line.geometry.setDrawRange(0, Math.floor(rb.nLine * lead)); rb.label.visible = lead >= .999; i++; }
    const grow = reduced ? 1 : Math.max(0, Math.min(1, (T - .7) / .25)); barMeshes.forEach((b, k) => { const g = Math.max(.001, Math.min(1, grow * 1.3 - k * .1)); b.m.scale.y = g; b.m.position.y = Y(3.5) + b.h * g / 2; b.l.visible = g >= .999; }); barTitle.visible = grow > .2;
    st.setBloom(.3 + .25 * S.hero.glow); if (st.bokeh) st.bokeh.uniforms.focus.value = camera.position.distanceTo(ctr.target);
    if (Math.abs(T - lastT) > 1e-4) { lastT = T; changed = true; }
    return changed;
  });
})();
