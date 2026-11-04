/* Body scenes for the LMO page, drawn by the shared Viz renderer (assets/viz.js) from window.LMO (assets/data.js):
   1. cell: the coin-cell stack exploded (Li disc, separator, electrolyte, LiMn2O4 coating, Al foil) to the measured thicknesses;
   2. curves: the four measured discharge curves as ribbons with the DFN model threaded through each;
   3. cv: the three cyclic voltammograms, measured and modelled, as ribbons over potential;
   4. loss: the DFN voltage-loss components along the 2C discharge as stacked ribbons. */
import { Viz } from './viz.js';
import { THREE } from './three-kit.js';
const D = window.LMO; const f = (v, d = 0) => (+v).toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d });
const RATE = { 'C/10': 'cyan', 'C/5': 'green', 'C/2': 'warm', '2C': 'rose' };

/* ---------- 1. the half cell, exploded ---------- */
Viz.view('[data-viz="cell"]', v => {
  const C = D.cell || { coating_um: 85, foil_um: 15, particle_um: 5, disc_mm: 10 }; const R = 1.5; const disc = (h, key, o = {}) => new THREE.Mesh(new THREE.CylinderGeometry(R, R, h, 64), Viz.mat(key, o));
  const layers = [['Li metal disc', .34, 'steel', { metalness: .8, roughness: .3 }], ['separator, electrolyte-filled', .12, 'white', { opacity: .55 }], [`LiMn₂O₄ coating, ${C.coating_um} µm`, .3 * C.coating_um / 85, 'violet', { roughness: .6 }], [`Al foil, ${C.foil_um} µm`, .06, 'steel', { metalness: .9, roughness: .25 }]];
  let y = -.5; const parts = [], labs = [];
  layers.forEach(([name, h, key, o], i) => { const m = disc(h, key, o); m.position.y = y + h / 2; m.userData = { name }; v.pickables.push(m); parts.push({ mesh: m, from: [0, 1.8 + i * .5, 0], delay: i * .16, span: .3, fade: false }); const ly = -.75 + i * .42; labs.push(Viz.label(v, name, [R + 1.25, ly, 0], { px: 20, scale: .09, color: 'ink' })); const ld = Viz.line(v, [[R * .95, y + h / 2, 0], [R + .55, ly, 0]], 'muted'); ld.material.transparent = true; ld.material.opacity = .7; labs[labs.length - 1].userData.lead = ld; y += h; });
  // particles in the coating: a ring of spheres at the measured particle size ratio
  const pg = new THREE.SphereGeometry(.045, 10, 8), pm = Viz.mat('violet', { roughness: .5 }); const NP = 160; const inst = new THREE.InstancedMesh(pg, pm, NP); const m4 = new THREE.Matrix4(); let s = 3; const rnd = () => { s = (s * 16807) % 2147483647; return (s - 1) / 2147483646; };
  const cy = parts[2].mesh.position.y, ch = layers[2][1]; for (let i = 0; i < NP; i++) { const r = Math.sqrt(rnd()) * (R - .06), th = rnd() * Math.PI * 2; m4.makeScale(1, 1, 1).setPosition(Math.cos(th) * r, cy + ch / 2 + .03, Math.sin(th) * r); inst.setMatrixAt(i, m4); } parts[2].mesh.add(inst); inst.position.y = -cy;
  const ex = Viz.exploded(v, parts);
  const foot = Viz.label(v, `${C.disc_mm} mm disc · ${C.particle_um} µm particles · 3.5 to 4.5 V`, [0, -.95, R + .2], { px: 20, scale: .09, bg: false });
  v.update = vv => { ex.draw(vv.p); labs.forEach((l, i) => { const o = Math.max(0, Math.min(1, (vv.p - i * .16 - .2) * 4)); l.material.opacity = o; l.userData.lead.material.opacity = o * .7; }); foot.material.opacity = Math.max(0, vv.p * 2 - 1); return false; };
  v.hover = hit => hit.object.userData && hit.object.userData.name ? `<b>${hit.object.userData.name}</b>` : null;
  v.spin = 0; v.rot = -.5; v.tilt = .35;
}, { camera: [.6, 2.4, 6.0], target: [.5, -.1, 0], fov: 34 });

/* ---------- 2. discharge curves ---------- */
Viz.view('[data-viz="curves"]', v => {
  const A = Viz.axes(v, { size: [4.6, 2.2, 2.6], x: { lo: 0, hi: 1.5, ticks: [0, .5, 1, 1.5], label: 'discharge capacity, mAh' }, y: { lo: 3.4, hi: 4.6, ticks: [3.5, 4, 4.5], label: 'V' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const rbs = D.curves.map((c, i) => { const rb = Viz.ribbon(v, A, c.measured.map(p => p[0]), c.measured.map(p => p[1]), { color: RATE[c.rate], width: .18, z: i, name: `${c.rate} measured`, unit: 'V' }); rb.mesh.userData.curve = c; const ln = Viz.line(v, c.model.map(p => { const q = A.map(p[0], p[1] + .006, i); return [q.x, q.y, q.z]; }), 'ink'); ln.material.transparent = true; ln.material.opacity = .9; ln.geometry.setDrawRange(0, 0); ln.userData.n = c.model.length; const lab = Viz.label(v, `${c.rate}${c.role === 'blind prediction' ? ' · blind' : c.role === 'training' ? ' · fit' : ''}`, [A.map(1.52, c.measured[c.measured.length - 1][1], i).x + .35, A.map(0, c.measured[c.measured.length - 1][1]).y, A.map(0, 0, i).z], { px: 20, scale: .09, color: RATE[c.rate], weight: 700 }); return { rb, ln, lab }; });
  v.update = vv => { rbs.forEach((r, i) => { const p = Math.max(0, Math.min(1, vv.p * 1.3 - i * .08)); r.rb.draw(p); r.ln.geometry.setDrawRange(0, Math.floor(r.ln.userData.n * p)); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const c = hit.object.userData.curve; if (!c) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${c.rate}, ${c.role}</b>measured ${f(c.measured[i][1], 3)} V at ${f(c.measured[i][0], 3)} mAh · Q ${c.q_meas} / ${c.q_model} mAh · V-Q RMSE ${c.vq_rmse_mV} mV`; };
  v.rot = -.4; v.tilt = .12;
}, { camera: [0, 2.2, 7.6], target: [0, -.2, 0], fov: 32 });

/* ---------- 3. cyclic voltammetry ---------- */
Viz.view('[data-viz="cv"]', v => {
  const A = Viz.axes(v, { size: [4.6, 2.2, 2.2], x: { lo: 3.0, hi: 4.9, ticks: [3.5, 4, 4.5], label: 'potential, V' }, y: { lo: -3.5, hi: 3.8, ticks: [-3, 0, 3], label: 'mA' }, z: { lo: -.5, hi: 2.5, ticks: [] } });
  const cols = ['cyan', 'green', 'warm'];
  const sets = D.cv_curves.map((c, i) => { const rb = Viz.ribbon(v, A, c.measured.map(p => p[0]), c.measured.map(p => p[1]), { color: cols[i], width: .14, z: i, name: `${c.scan} measured`, skirt: false }); rb.mesh.userData.cv = c; const ln = Viz.line(v, c.model.map(p => { const q = A.map(p[0], p[1], i); return [q.x, q.y + .01, q.z]; }), 'ink'); ln.material.transparent = true; ln.material.opacity = .85; ln.geometry.setDrawRange(0, 0); ln.userData.n = c.model.length; const lab = Viz.label(v, c.scan, [A.map(4.9, 0, i).x + .4, A.map(0, 3.0).y, A.map(0, 0, i).z], { px: 20, scale: .09, color: cols[i], weight: 700 }); return { rb, ln, lab }; });
  Viz.label(v, 'solid: measured · line: DFN with the measured voltage programme', [0, A.map(0, 3.8).y + .3, -1.2], { px: 20, scale: .085, bg: false });
  v.update = vv => { sets.forEach((r, i) => { const p = Math.max(0, Math.min(1, vv.p * 1.3 - i * .1)); r.rb.draw(p); r.ln.geometry.setDrawRange(0, Math.floor(r.ln.userData.n * p)); r.lab.material.opacity = p > .95 ? 1 : 0; }); return false; };
  v.hover = hit => { const c = hit.object.userData.cv; if (!c) return null; const i = Viz.ribbonIndex(hit, A); const row = D.cv_rows.find(r => r.scan === c.scan); return `<b>${c.scan}</b>${f(c.measured[i][1], 2)} mA at ${f(c.measured[i][0], 2)} V · peaks meas / DFN: ox ${row.ox_meas} / ${row.ox_model} V, red ${row.red_meas} / ${row.red_model} V`; };
  v.rot = -.35; v.tilt = .1;
}, { camera: [0, 2.2, 7.6], target: [0, -.2, 0], fov: 32 });

/* ---------- 4. voltage-loss components along the 2C discharge ---------- */
const lossView = Viz.view('[data-viz="loss"]', v => {
  const P = D.op_profile; const q = P.capacity_mAh; const comps = [['contact_ohmic', 'contact / ohmic', 'rose'], ['reaction', 'charge transfer', 'warm'], ['particle_diffusion', 'particle diffusion', 'green'], ['electrolyte', 'electrolyte', 'cyan']];
  const ely = P.electrolyte_concentration.map((x, i) => x + P.electrolyte_ohmic[i] + P.solid_ohmic[i]); const series = { contact_ohmic: P.contact_ohmic, reaction: P.reaction, particle_diffusion: P.particle_diffusion.map(Math.abs), electrolyte: ely.map(Math.abs) };
  const ymax = Math.max(...q.map((_, i) => comps.reduce((a, c) => a + series[c[0]][i], 0))) * 1.1;
  const A = Viz.axes(v, { size: [4.4, 2.2, 2.2], x: { lo: 0, hi: Math.max(...q), ticks: [0, .25, .5, .75], fmt: t => t.toFixed(2), label: 'capacity, mAh' }, y: { lo: 0, hi: ymax, ticks: [0, 100, 200, 300, 400], label: 'mV' }, z: { lo: -.5, hi: 3.5, ticks: [] } });
  const rbs = comps.map(([k, name, col], i) => { const ys = series[k]; const rb = Viz.ribbon(v, A, q, ys, { color: col, width: .34, z: i, name, unit: 'mV' }); rb.mesh.userData.k = k; const lab = Viz.label(v, `${name} · ${f(D.stats.shares[k] ?? (D.stats.shares.electrolyte_concentration + D.stats.shares.electrolyte_ohmic + D.stats.shares.solid_ohmic), 0)} %`, [A.map(Math.max(...q), 0, i).x + .5, A.map(0, ys[ys.length - 1]).y + .1, A.map(0, 0, i).z], { px: 20, scale: .09, color: col, weight: 700 }); return { rb, lab, i }; });
  v.step = 0;
  v.update = vv => { rbs.forEach((r, i) => { const on = vv.step === 0 || vv.step > 3 || (vv.step === 1 && i === 0) || (vv.step === 2 && i === 1) || (vv.step === 3 && i >= 2); r.rb.mesh.material.opacity = on ? 1 : .25; r.rb.mesh.material.transparent = true; r.rb.draw(vv.p); r.lab.material.opacity = vv.p > .9 ? (on ? 1 : .3) : 0; }); return false; };
  v.hover = hit => { const k = hit.object.userData.k; if (!k) return null; const i = Viz.ribbonIndex(hit, A); return `<b>${hit.object.userData.series} at ${f(q[i], 3)} mAh</b>${f(series[k][i])} mV · terminal ${f(P.voltage_V[i], 3)} V, OCV ${f(P.ocv_V[i], 3)} V`; };
  v.rot = -.4; v.tilt = .1;
}, { camera: [0, 2.0, 7.2], target: [0, -.2, 0], fov: 32 });
document.addEventListener('lmo-step', e => { if (lossView) { lossView.step = +e.detail; lossView.dirty = true; } });
