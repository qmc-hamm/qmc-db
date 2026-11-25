import { fmtErr, fmtET, keepAll } from "./fmt.js";

const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
const PAGE = 100;
const SCAT = (() => { try { return document.createElement("canvas").getContext("webgl") ? "scattergl" : "scatter"; } catch (e) { return "scatter"; } })();  // scattergl needs WebGL
const PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f", "#d62728"];
const TABS = [["hydrogen", "Hydrogen"], ["graphene", "Graphene"], ["defects", "Defects"], ["li7p3s11", "Li7P3S11"]];
const FILTERS = [["source", "source"], ["method", "method"], ["creation_path", "creation_path"], ["phase", "phase"], ["natoms", "N"], ["qmc_quality", "qmc_quality"]];

let C = null;           // columnar index
let N = 0;
let meta = {};
let tab = "hydrogen";
let selection = null;   // Set of row indices from plot box/lasso, or null
let noTP = false;       // table restricted to rows lacking T or P_used
let filtered = [];      // row indices passing filters (current tab)
let view = [];          // filtered ∩ selection (∩ noTP)
let sortKey = null, sortDir = 1, page = 0;
let curUuid = null;
let uuidMap = null;

const rs = (i) => C.rs_computed[i] ?? C.rs_recorded[i];
const eatom = (i) => (C.E_per_atom[i] ?? (C.E[i] != null && C.natoms[i] ? C.E[i] / C.natoms[i] : null));
const eatomErr = (i) => (C.E_err[i] != null && C.natoms[i] ? C.E_err[i] / C.natoms[i] : null);
const hasTP = (i) => C.T[i] != null && C.P_used[i] != null;

const HCOLS = [
  ["uuid", (i) => `<span title="${esc(C.uuid[i])}">${esc(C.uuid[i].slice(0, 8))}</span>`, (i) => C.uuid[i]],
  ["source", (i) => esc(C.source[i]), (i) => C.source[i]],
  ["method", (i) => esc(C.method[i]), (i) => C.method[i]],
  ["N", (i) => esc(C.natoms[i]), (i) => C.natoms[i]],
  ["T (K)", (i) => esc(C.T[i] ?? "—"), (i) => C.T[i]],
  ["P (GPa)", (i) => (C.P_used[i] == null ? "—" : C.P_used[i].toFixed(1) + (C.P_extrapolated[i] ? " (extrap.)" : C.P_estimated[i] ? " (est.)" : "")), (i) => C.P_used[i]],
  ["r_s", (i) => (rs(i) == null ? "—" : rs(i).toFixed(3)), (i) => rs(i)],
  ["E/atom (eV)", (i) => fmtErr(eatom(i), eatomErr(i)), (i) => eatom(i)],
  ["phase", (i) => esc(C.phase[i]), (i) => C.phase[i]],
  ["⬇", (i) => `<a href="/h5/${encodeURIComponent(C.uuid[i])}.h5" title="download .h5">⬇</a>`, null],
];
const OCOLS = [
  ["system", (i) => esc(C.system[i]), (i) => C.system[i]],
  ["formula", (i) => esc(C.formula[i]), (i) => C.formula[i]],
  ["method", (i) => esc(C.method[i]), (i) => C.method[i]],
  ["N", (i) => esc(C.natoms[i]), (i) => C.natoms[i]],
  ["E (eV)", (i) => fmtErr(C.E[i], C.E_err[i]), (i) => C.E[i]],
  ["author", (i) => esc([].concat(C.author[i] ?? "").join(", ")), (i) => String([].concat(C.author[i] ?? ""))],
  ["⬇", (i) => `<a href="/h5/${encodeURIComponent(C.uuid[i])}.h5" title="download .h5">⬇</a>`, null],
];

async function init() {
  const r = await fetch("/data/index.json.gz");
  const idx = await r.json();
  C = idx.columns; N = C.uuid.length; meta = idx;
  header();
  buildTabs();
  buildFilters();
  $("colorby").onchange = $("xaxis").onchange = () => { selection = null; page = 0; refresh(true); };
  $("clearsel").onclick = () => { selection = null; refresh(true); };
  $("dl").onclick = download;
  $("notp").onclick = () => { noTP = !noTP; page = 0; refresh(false); };
  $("eos").addEventListener("toggle", loadEos, { once: true });
  setTab("hydrogen");
}

function header() {
  const bySrc = {};
  for (let i = 0; i < N; i++) bySrc[C.source[i] ?? "none"] = (bySrc[C.source[i] ?? "none"] || 0) + 1;
  const parts = Object.entries(bySrc).sort((a, b) => b[1] - a[1]).map(([k, v]) => `${esc(k)}: ${v}`);
  $("hdr").innerHTML = `${N} files · ${parts.join(", ")} · harvested ${fmtET(meta.harvested_at)} · EOS: ${esc(meta.eos)} (backup: none)`;
}

function buildTabs() {
  $("tabs").innerHTML = "";
  for (const [k, label] of TABS) {
    const b = document.createElement("button");
    b.textContent = label; b.dataset.tab = k;
    b.onclick = () => setTab(k);
    $("tabs").append(b);
  }
}

function setTab(k) {
  tab = k; selection = null; noTP = false; page = 0; sortKey = null; curUuid = null;
  for (const b of $("tabs").children) b.classList.toggle("on", b.dataset.tab === k);
  $("tab-hydrogen").hidden = k !== "hydrogen";
  $("tab-other").hidden = k === "hydrogen";
  $("detail").hidden = true;
  $("othernote").textContent = k === "li7p3s11" ? "Li7P3S11 files are indexed from the repo, not the bucket." : "";
  refresh(true);
}

function buildFilters() {
  const bar = $("filters"); bar.innerHTML = "";
  for (const [f, label] of FILTERS) {
    const vals = new Map();
    for (let i = 0; i < N; i++) if (C.tab[i] === "hydrogen") { const v = String(C[f][i] ?? "∅"); vals.set(v, (vals.get(v) || 0) + 1); }
    const opts = [...vals.keys()].sort((a, b) => (isNaN(a) || isNaN(b) ? a.localeCompare(b) : a - b));
    const l = document.createElement("label");
    l.innerHTML = `${label}<select multiple id="f_${f}"><option class="all" selected>all</option>${opts.map((v) => `<option value="${esc(v)}">${esc(v)} (${vals.get(v)})</option>`).join("")}</select>`;
    bar.append(l);
    const sel = l.firstElementChild, all = sel.options[0];
    let wasAll = true;
    sel.onchange = () => {
      const others = [...sel.options].slice(1);
      if (keepAll(wasAll, all.selected, others.some((o) => o.selected))) { others.forEach((o) => { o.selected = false; }); all.selected = true; }
      else all.selected = false;
      wasAll = all.selected;
      selection = null; page = 0; refresh(true);
    };
  }
  for (const [id, label] of [["Tmin", "T min"], ["Tmax", "T max"], ["Pmin", "P min"], ["Pmax", "P max"]]) {
    const l = document.createElement("label");
    l.innerHTML = `${label}<input type="number" id="${id}">`;
    bar.append(l);
    l.lastElementChild.onchange = () => { selection = null; page = 0; refresh(true); };
  }
}

function computeFiltered() {
  const out = [];
  if (tab !== "hydrogen") {
    for (let i = 0; i < N; i++) if (C.tab[i] === tab) out.push(i);
    return out;
  }
  const sets = FILTERS.map(([f]) => {
    const s = [...$("f_" + f).selectedOptions].filter((o) => o.index > 0).map((o) => o.value);
    return s.length ? [C[f], new Set(s)] : null;
  }).filter(Boolean);
  const num = (id) => { const v = $(id).value; return v === "" ? null : Number(v); };
  const [tmin, tmax, pmin, pmax] = ["Tmin", "Tmax", "Pmin", "Pmax"].map(num);
  for (let i = 0; i < N; i++) {
    if (C.tab[i] !== "hydrogen") continue;
    let ok = true;
    for (const [col, s] of sets) if (!s.has(String(col[i] ?? "∅"))) { ok = false; break; }
    if (!ok) continue;
    const T = C.T[i], P = C.P_used[i];
    if ((tmin != null && (T == null || !(T >= tmin))) || (tmax != null && (T == null || !(T <= tmax))) || (pmin != null && (P == null || !(P >= pmin))) || (pmax != null && (P == null || !(P <= pmax)))) continue;
    out.push(i);
  }
  return out;
}

function refresh(redrawPlots) {
  filtered = computeFiltered();
  if (redrawPlots && tab === "hydrogen") drawPlots();
  updateView();
}

function updateView() {
  let v = filtered;
  if (tab === "hydrogen") {
    if (selection) v = v.filter((i) => selection.has(i));
    if (noTP) v = v.filter((i) => !hasTP(i));
    notpBar();
    $("selinfo").textContent = selection ? `${selection.size} selected` : "no selection";
    $("clearsel").disabled = !selection;
  }
  view = v;
  renderTable();
}

function notpBar() {
  const miss = filtered.filter((i) => !hasTP(i));
  const srcs = {};
  for (const i of miss) srcs[C.source[i] ?? "none"] = (srcs[C.source[i] ?? "none"] || 0) + 1;
  const b = $("notp");
  b.classList.toggle("on", noTP);
  b.textContent = `No T/P recorded: ${miss.length} files (${Object.entries(srcs).map(([k, v]) => `${k}: ${v}`).join(", ") || "none"})` + (noTP ? " — table filtered, click to undo" : "");
}

// ---------- plots ----------
function themeLayout(extra) {
  const cs = getComputedStyle(document.documentElement);
  const fg = cs.getPropertyValue("--fg").trim(), line = cs.getPropertyValue("--line").trim();
  return Object.assign({
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)", font: { color: fg },
    margin: { t: 20, r: 10, b: 50, l: 60 }, dragmode: "select",
    xaxis: { gridcolor: line, zerolinecolor: line }, yaxis: { gridcolor: line, zerolinecolor: line },
  }, extra);
}

function drawPlots() { drawPT(); drawERs(); }

function drawPT() {
  const useRs = $("xaxis").value === "rs";
  const by = $("colorby").value;
  const colourOf = new Map();
  const groups = new Map();
  const add = (key, sym, colour, i, hov) => {
    let g = groups.get(key);
    if (!g) groups.set(key, g = { x: [], y: [], cd: [], text: [], sym, colour, ext: key.endsWith("(P extrapolated)") });
    g.x.push(useRs ? rs(i) : C.P_used[i]); g.y.push(C.T[i]); g.cd.push(i); g.text.push(hov);
  };
  for (const i of filtered) {
    const x = useRs ? rs(i) : C.P_used[i];
    if (x == null || C.T[i] == null) continue;
    const gname = String(C[by][i] ?? "none");
    if (!colourOf.has(gname)) colourOf.set(gname, PALETTE[colourOf.size % PALETTE.length]);
    const col = colourOf.get(gname);
    let hov = `${esc(C.uuid[i].slice(0, 8))} ${esc(C.source[i])} ${esc(C.method[i])} N=${esc(C.natoms[i])}<br>T=${esc(C.T[i])} K P=${C.P_used[i] == null ? "—" : C.P_used[i].toFixed(1)} GPa r_s=${rs(i)?.toFixed(3) ?? "—"}`;
    if (useRs && C.rs_suspect[i]) { add("r_s suspect", "circle-open", "#d62728", i, hov + "<br>r_s suspect"); continue; }
    if (!useRs && C.P_extrapolated[i]) { add(gname + " (P extrapolated)", "circle-open", col, i, hov + "<br>P estimated by EOS, extrapolated outside fit range"); continue; }
    if (!useRs && C.P_estimated[i]) { add(gname + " (P est.)", "circle-open", col, i, hov + "<br>P estimated from EOS"); continue; }
    add(gname, "circle", col, i, hov);
  }
  const traces = [...groups].map(([name, g]) => ({
    type: SCAT, mode: "markers", name, x: g.x, y: g.y, customdata: g.cd, text: g.text, hoverinfo: "text",
    marker: { symbol: g.sym, color: g.colour, size: g.sym === "circle" ? 6 : g.ext ? 12 : 8, line: { width: g.ext ? 3 : 1.5, color: g.colour } },
  }));
  const gd = $("pt");
  Plotly.react(gd, traces, themeLayout({
    xaxis: { title: useRs ? "r_s (bohr)" : "P (GPa)", gridcolor: "#8884" },
    yaxis: { title: "T (K)", gridcolor: "#8884" },
  }), { responsive: true });
  if (!gd._wired) { wire(gd); gd._wired = true; }
}

function drawERs() {
  const x = [], y = [], e = [], T = [], cd = [];
  for (const i of filtered) {
    if (C.method[i] !== "DMC" || rs(i) == null || eatom(i) == null) continue;
    x.push(rs(i)); y.push(eatom(i)); e.push(eatomErr(i) ?? 0); T.push(C.T[i]); cd.push(i);
  }
  const gd = $("ers");
  Plotly.react(gd, [{
    type: SCAT, mode: "markers", x, y, customdata: cd, error_y: { type: "data", array: e, visible: true, thickness: 1 },
    marker: { color: T, colorscale: "Viridis", size: 6, showscale: true, colorbar: { title: "T (K)" } },
    text: cd.map((i, k) => `${esc(C.uuid[i].slice(0, 8))} T=${esc(T[k])} E/atom=${fmtErr(y[k], e[k])}`), hoverinfo: "text",
  }], themeLayout({ xaxis: { title: "r_s (bohr)" }, yaxis: { title: cd.length && cd.every((i) => C.fsc_applied[i] === true) ? "E/atom (eV, FSC-corrected)" : "E/atom (eV)" } }), { responsive: true });
  if (!gd._wired) { wire(gd); gd._wired = true; }
}

function wire(gd) {
  gd.on("plotly_selected", (ev) => {
    if (!ev || !ev.points || !ev.points.length) return;
    selection = new Set(ev.points.map((p) => p.customdata));
    page = 0; updateView();
  });
  gd.on("plotly_deselect", () => { selection = null; page = 0; updateView(); });
}

// ---------- table ----------
function renderTable() {
  const cols = tab === "hydrogen" ? HCOLS : OCOLS;
  if (sortKey != null) {
    const f = cols[sortKey][2];
    view = [...view].sort((a, b) => {
      const x = f(a), y = f(b);
      if (x == null) return y == null ? 0 : 1;   // nulls last
      if (y == null) return -1;
      return (x < y ? -1 : x > y ? 1 : 0) * sortDir;
    });
  }
  const np = Math.max(1, Math.ceil(view.length / PAGE));
  page = Math.min(page, np - 1);
  const rows = view.slice(page * PAGE, (page + 1) * PAGE);
  const head = cols.map((c, k) => `<th data-k="${k}">${c[0]}${sortKey === k ? (sortDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("");
  const body = rows.map((i) => `<tr data-i="${i}" class="${C.uuid[i] === curUuid ? "cur" : ""}">${cols.map((c) => `<td>${c[1](i)}</td>`).join("")}</tr>`).join("");
  $("tbl").innerHTML = `<thead><tr>${head}</tr></thead><tbody>${body}</tbody>`;
  for (const th of $("tbl").querySelectorAll("th")) th.onclick = () => {
    const k = +th.dataset.k; if (!cols[k][2]) return;
    sortDir = sortKey === k ? -sortDir : 1; sortKey = k; renderTable();
  };
  for (const tr of $("tbl").querySelectorAll("tbody tr")) tr.onclick = (ev) => { if (ev.target.tagName !== "A") showDetail(C.uuid[+tr.dataset.i]); };
  $("tablebar").innerHTML = `<b>${view.length}</b> rows` + (tab === "hydrogen" ? " (filters ∩ selection)" : "");
  $("pager").innerHTML = `<button id="pp" ${page === 0 ? "disabled" : ""}>prev</button> page ${page + 1} / ${np} <button id="pn" ${page >= np - 1 ? "disabled" : ""}>next</button>`;
  $("pp").onclick = () => { page--; renderTable(); };
  $("pn").onclick = () => { page++; renderTable(); };
}

// ---------- detail ----------
// S(k): either structure_factor_k (Nk,) or |k| from structure_factor_ks (Nk,3); S may be (Nk,) or (Nk,3).
function skData(ob) {
  const S = ob.structure_factor;
  if (!Array.isArray(S) || !S.length) return null;
  const col = (a) => (Array.isArray(a) ? a.map((v) => (Array.isArray(v) ? v[0] : v)) : null);
  let k = Array.isArray(ob.structure_factor_k) && ob.structure_factor_k.length ? ob.structure_factor_k : null;
  if (!k && Array.isArray(ob.structure_factor_ks) && ob.structure_factor_ks.length)
    k = ob.structure_factor_ks.map((v) => (Array.isArray(v) ? Math.hypot(...v) : Math.abs(v)));
  if (!k || k.length !== S.length) return null;
  const twoD = Array.isArray(S[0]);
  const err = col(ob.structure_factor_error);
  return { k, y: col(S), err: err && err.length === S.length ? err : null,
    xlabel: "k (1/bohr)", ylabel: twoD ? "S(k) (first column)" : "S(k)" };
}

function getUuidMap() {
  if (uuidMap == null) { uuidMap = new Map(); for (let i = 0; i < N; i++) uuidMap.set(C.uuid[i], i); }
  return uuidMap;
}

async function showDetail(uuid) {
  curUuid = uuid;
  const box = $("detail"); box.hidden = false; box.textContent = "loading " + uuid + "…";
  let d;
  try {
    const r = await fetch("/api/file/" + encodeURIComponent(uuid));
    if (!r.ok) throw new Error(r.status + " " + (await r.text()));
    d = await r.json();
  } catch (e) { box.innerHTML = `<span class="warn">${esc(e.message)}</span> (file not in local cache?)`; return; }
  if (curUuid !== uuid) return;
  const i = getUuidMap().get(uuid);
  const ob = d.observables || {};
  const lk = d.links || {};
  const lkl = (arr) => (arr || []).map((u) => `<span class="link" data-u="${esc(u)}">${esc(String(u).slice(0, 8))}</span>`).join(" ") || "—";
  const rowsOf = (o) => Object.entries(o).map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(typeof v === "object" ? JSON.stringify(v) : v)}</td></tr>`).join("");
  const fsc = Object.entries(ob).filter(([k]) => /fsc/i.test(k));
  const energy = `<tr><td>total energy (eV${i != null && C.fsc_applied[i] === true ? ", FSC-corrected" : ""})</td><td>${fmtErr(ob.total_energy, ob.total_energy_error)}</td></tr>`
    + `<tr><td>uncorrected (eV, not used as energy)</td><td>${fmtErr(ob.total_energy_uncorrected ?? (i != null ? C.E_unc[i] : null), null)}</td></tr>`
    + fsc.map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(typeof v === "number" ? v.toFixed(4) : JSON.stringify(v))}</td></tr>`).join("");
  const plain = Object.fromEntries(Object.entries(ob).filter(([k]) => !/^structure_factor/.test(k)));
  box.innerHTML = `<h3>${esc(uuid)} <a href="/h5/${encodeURIComponent(uuid)}.h5">⬇ .h5</a></h3>
    <p>config_id: ${esc(lk.config_id)} · linked DFT: ${lkl(lk.linked_dft_uuids)} · linked QMC: ${lkl(lk.linked_qmc_uuids)}
    ${(lk.missing_links || []).length ? `<span class="warn">· missing links: ${esc(lk.missing_links.length)}</span>` : ""}
    ${(lk.geom_mismatch || []).length ? `<span class="warn">· geometry mismatch</span>` : ""}</p>
    <h4>Energy</h4><table>${energy}</table>
    <h4>Attributes</h4><table>${rowsOf(d.attrs || {})}</table>
    <h4>Observables</h4><table>${rowsOf(plain)}</table>
    <h4>Parameters</h4><table>${rowsOf({ ...(d.parameters?.attrs || {}), ...(d.parameters?.values || {}) })}</table>
    <h4>Provenance</h4><table>${rowsOf({ ...(d.provenance?.attrs || {}), ...(d.provenance?.values || {}) })}</table>
    <div id="skplot" class="plot small" hidden></div>`;
  for (const s of box.querySelectorAll(".link")) s.onclick = () => showDetail(s.dataset.u);
  const sk = skData(ob);
  if (sk) {
    const el = $("skplot"); el.hidden = false;
    Plotly.newPlot(el, [{ type: "scatter", mode: "lines+markers", x: sk.k, y: sk.y, error_y: sk.err ? { type: "data", array: sk.err, visible: true } : undefined, marker: { size: 4 } }],
      themeLayout({ xaxis: { title: sk.xlabel }, yaxis: { title: sk.ylabel }, dragmode: "zoom" }), { responsive: true });
  }
  box.scrollIntoView({ behavior: "smooth", block: "nearest" });
  for (const tr of $("tbl").querySelectorAll("tbody tr")) tr.classList.toggle("cur", C.uuid[+tr.dataset.i] === uuid);
}

// ---------- download ----------
function exportPlan(rows) {
  let idx = [...rows];
  if ($("inclinked").checked) {
    const seen = new Set(idx), um = getUuidMap();
    for (const i of rows) for (const u of C.linked_dft_uuids[i] || []) {
      const j = um.get(u);
      if (j != null && !seen.has(j)) { seen.add(j); idx.push(j); }
    }
  }
  const inclSuspect = $("inclsuspect").checked;
  const suspect = idx.filter((i) => C.rs_suspect[i] === true).length;
  if (!inclSuspect) idx = idx.filter((i) => C.rs_suspect[i] !== true);
  const methods = {};
  for (const i of idx) methods[C.method[i]] = (methods[C.method[i]] || 0) + 1;
  const n = (f) => idx.filter(f).length;
  return { frames: idx.length, methods, suspect, inclSuspect,
    noDft: n((i) => C.method[i] !== "DFT" && !(C.linked_dft_uuids[i] || []).length),
    mism: n((i) => (C.geom_mismatch[i] || []).length), miss: n((i) => (C.missing_links[i] || []).length),
    eNull: n((i) => C.E[i] == null) };
}

function download() {
  const rows = view;
  if (!rows.length) { alert("Nothing to download."); return; }
  const p = exportPlan(rows);
  const ms = Object.entries(p.methods);
  if (ms.length > 1 || p.noDft || p.mism || p.miss || p.suspect || p.eNull) {
    const msg = `Frames planned: ${p.frames}\nMethods: ${ms.map(([k, v]) => `${k} ×${v}`).join(", ") || "none"}\n`
      + `QMC rows with no linked DFT: ${p.noDft}\nRows with geometry mismatch: ${p.mism}\nRows with missing links: ${p.miss}\n`
      + `rs_suspect rows: ${p.suspect} (${p.inclSuspect ? "included" : "skipped"})\nRows with energy missing (NaN): ${p.eNull}\nExport?`;
    if (!confirm(msg)) return;
  }
  if (!selection && rows.length > 5000 && !confirm(`No selection: export all ${rows.length} filtered rows?`)) return;
  $("xpayload").value = JSON.stringify({ uuids: rows.map((i) => C.uuid[i]), include_linked_dft: $("inclinked").checked, include_suspect: $("inclsuspect").checked });
  $("xform").submit();
}

// ---------- EOS ----------
async function loadEos() {
  const el = $("eosbody");
  let recs;
  try {
    const r = await fetch("/data/eos_residuals.json");
    if (!r.ok) throw new Error(r.status);
    recs = await r.json();
  } catch (e) { el.textContent = "No EOS residuals available (eos: " + meta.eos + ")."; return; }
  if (!recs.length) { el.textContent = "EOS residual list is empty."; return; }
  const med = (a) => { const s = [...a].sort((x, y) => x - y), n = s.length; return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2; };
  recs = recs.filter((r) => r.rs != null && r.T != null && r.rel_resid != null);
  if (!recs.length) { el.textContent = "EOS residual list has no usable rows."; return; }
  let lo = Infinity, hi = -Infinity;
  for (const r of recs) { if (r.rs < lo) lo = r.rs; if (r.rs > hi) hi = r.rs; }
  const nb = 20, w = (hi - lo) / nb || 1;
  const Ts = [...new Set(recs.map((r) => r.T))].sort((a, b) => a - b);
  const cells = new Map();
  for (const r of recs) {
    const b = Math.min(nb - 1, Math.floor((r.rs - lo) / w)), key = b + "|" + r.T;
    (cells.get(key) || cells.set(key, []).get(key)).push(r.rel_resid);
  }
  const z = Ts.map((T) => Array.from({ length: nb }, (_, b) => (cells.has(b + "|" + T) ? med(cells.get(b + "|" + T)) : null)));
  const groups = {};
  for (const r of recs) (groups[r.source + " / " + r.creation_path] ||= []).push(r.rel_resid);
  const trs = Object.entries(groups).map(([k, a]) => `<tr><td>${esc(k)}</td><td>${a.length}</td><td>${(100 * med(a)).toFixed(2)}%</td><td>${(100 * a.reduce((m, x) => Math.max(m, Math.abs(x)), 0)).toFixed(2)}%</td></tr>`).join("");
  el.innerHTML = `<div id="eosmap" class="plot small"></div><table><thead><tr><th>source / creation_path</th><th>n</th><th>median rel. resid.</th><th>max |rel. resid.|</th></tr></thead><tbody>${trs}</tbody></table>`;
  Plotly.newPlot("eosmap", [{ type: "heatmap", z, x: Array.from({ length: nb }, (_, b) => +(lo + (b + 0.5) * w).toFixed(3)), y: Ts, colorscale: "RdBu", zmid: 0, colorbar: { title: "median rel_resid" }, hoverongaps: false }],
    themeLayout({ xaxis: { title: "r_s (bohr)" }, yaxis: { title: "T (K)", type: "category" }, dragmode: "zoom" }), { responsive: true });
}

init().catch((e) => { $("hdr").innerHTML = `<span class="warn">failed to load index: ${esc(e.message)}</span>`; });
