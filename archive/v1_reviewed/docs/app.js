(function () {
  "use strict";

  const DATA = window.MIZER_MODEL_DATA;
  if (!DATA) {
    console.error("MIZER_MODEL_DATA not found. Run R/03_generate_figures_and_web_data.R first.");
    return;
  }

  // State
  const state = {
    fp: 0.15,
    fm: 0.05,
    fcod: 0.00,
    taper: false,
    activePreset: "status_quo",
    spectrumView: "spectra",      // "spectra" | "mortality" | "feeding"
    regimeCompare: "moratorium"   // "moratorium" | "cultivated" | "both"
  };

  const PRESETS = DATA.presets || {
    status_quo:        { fp: 0.15, fm: 0.05, fcod: 0.00, taper: false },
    pelagic_rebalance: { fp: 0.55, fm: 0.05, fcod: 0.00, taper: false },
    ecosystem_pulse:   { fp: 0.65, fm: 0.45, fcod: 0.15, taper: true  },
    premature_reopen:  { fp: 0.15, fm: 0.05, fcod: 0.25, taper: false }
  };

  // Populate Hero Stats from exact Mizer model metrics
  if (DATA.hero) {
    const setTxt = (id, val) => {
      const el = document.getElementById(id);
      if (el && val !== undefined) el.innerHTML = val;
    };
    setTxt("hero-collapse-pct", DATA.hero.collapse_pct);
    setTxt("hero-forage-surge", DATA.hero.forage_surge);
    setTxt("hero-larval-mort", DATA.hero.larval_mort_ratio);
    setTxt("hero-recovery-time", DATA.hero.recovery_window);
  }

  // Bilinear + trilinear interpolation over the exact precomputed Mizer scenario grid
  function getInterpolatedScenario(fp, fm, fcod, taper) {
    const grid = taper ? DATA.grid_taper : DATA.grid_sustained;
    const fpVals = DATA.grid_axes.fp;
    const fmVals = DATA.grid_axes.fm;
    const fcVals = DATA.grid_axes.fcod;

    function bracket(arr, x) {
      if (x <= arr[0]) return [0, 0, 0];
      const n = arr.length - 1;
      if (x >= arr[n]) return [n, n, 0];
      for (let i = 0; i < n; i++) {
        if (x >= arr[i] && x <= arr[i + 1]) {
          const t = (x - arr[i]) / (arr[i + 1] - arr[i] || 1);
          return [i, i + 1, t];
        }
      }
      return [0, 0, 0];
    }

    const [i0, i1, tx] = bracket(fpVals, fp);
    const [j0, j1, ty] = bracket(fmVals, fm);
    const [k0, k1, tz] = bracket(fcVals, fcod);

    function key(i, j, k) {
      return `${i}_${j}_${k}`;
    }

    const c000 = grid[key(i0, j0, k0)];
    const c100 = grid[key(i1, j0, k0)];
    const c010 = grid[key(i0, j1, k0)];
    const c110 = grid[key(i1, j1, k0)];
    const c001 = grid[key(i0, j0, k1)];
    const c101 = grid[key(i1, j0, k1)];
    const c011 = grid[key(i0, j1, k1)];
    const c111 = grid[key(i1, j1, k1)];

    function lerpArr(field) {
      const len = c000[field].length;
      const out = new Array(len);
      for (let m = 0; m < len; m++) {
        const v00 = c000[field][m] * (1 - tx) + c100[field][m] * tx;
        const v10 = c010[field][m] * (1 - tx) + c110[field][m] * tx;
        const v01 = c001[field][m] * (1 - tx) + c101[field][m] * tx;
        const v11 = c011[field][m] * (1 - tx) + c111[field][m] * tx;
        const v0 = v00 * (1 - ty) + v10 * ty;
        const v1 = v01 * (1 - ty) + v11 * ty;
        out[m] = v0 * (1 - tz) + v1 * tz;
      }
      return out;
    }

    function lerpScalar(field) {
      const v00 = c000[field] * (1 - tx) + c100[field] * tx;
      const v10 = c010[field] * (1 - tx) + c110[field] * tx;
      const v01 = c001[field] * (1 - tx) + c101[field] * tx;
      const v11 = c011[field] * (1 - tx) + c111[field] * tx;
      const v0 = v00 * (1 - ty) + v10 * ty;
      const v1 = v01 * (1 - ty) + v11 * ty;
      return v0 * (1 - tz) + v1 * tz;
    }

    return {
      cod_ssb: lerpArr("cod_ssb"),
      forage_bio: lerpArr("forage_bio"),
      sculpin_bio: lerpArr("sculpin_bio"),
      cod_yield: lerpArr("cod_yield"),
      pelagic_yield: lerpArr("pelagic_yield"),
      other_yield: lerpArr("other_yield"),
      larval_mort_end: lerpScalar("larval_mort_end")
    };
  }

  // HiDPI Canvas helper
  function setupCanvas(canvas) {
    const dpr = window.devicePixelRatio || 1;
    const rect = canvas.getBoundingClientRect();
    const width = rect.width || canvas.width;
    const height = parseFloat(canvas.getAttribute("height")) || 300;
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, width, height };
  }

  // Draw Time Series Chart (Historical + Policy Projection)
  function drawTimeSeries(scen) {
    const canvas = document.getElementById("canvas-timeseries");
    if (!canvas) return;
    const { ctx, width, height } = setupCanvas(canvas);
    ctx.clearRect(0, 0, width, height);

    const hist = DATA.history;
    const projYears = DATA.proj_years; // 0..40 (mapped to calendar 2026..2066)
    const years = hist.years.concat(projYears.slice(1));
    const codSeries = hist.cod_ssb.concat(scen.cod_ssb.slice(1));
    const forageSeries = hist.forage_bio.concat(scen.forage_bio.slice(1));
    const sculpinSeries = hist.sculpin_bio.concat(scen.sculpin_bio.slice(1));

    const pad = { left: 56, right: 20, top: 22, bottom: 34 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const xMin = years[0];
    const xMax = years[years.length - 1];
    const rawMax = Math.max(...codSeries, ...forageSeries, DATA.reference.precollapse_cod_ssb * 1.15);
    const yMax = Math.ceil(rawMax / 200) * 200;

    const xScale = (yr) => pad.left + ((yr - xMin) / (xMax - xMin)) * pw;
    const yScale = (v) => pad.top + ph - (v / yMax) * ph;

    // Shaded historical phases
    const x1992 = xScale(1992);
    const xPolicy = xScale(DATA.policy_start_year);

    // Moratorium trap shading (1992 -> Policy Start)
    ctx.fillStyle = "rgba(197, 57, 41, 0.055)";
    ctx.fillRect(x1992, pad.top, xPolicy - x1992, ph);

    // Policy intervention shading (Policy Start -> end)
    ctx.fillStyle = "rgba(29, 95, 167, 0.05)";
    ctx.fillRect(xPolicy, pad.top, pad.left + pw - xPolicy, ph);

    // Grid lines
    ctx.strokeStyle = "#e5ecf3";
    ctx.lineWidth = 1;
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#5c7089";
    const yTicks = 5;
    for (let i = 0; i <= yTicks; i++) {
      const val = (yMax / yTicks) * i;
      const y = yScale(val);
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(pad.left + pw, y);
      ctx.stroke();
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(Math.round(val) + " kt", pad.left - 7, y);
    }

    // Cultivation threshold dashed line
    if (DATA.reference.allee_threshold_ssb) {
      const yThr = yScale(DATA.reference.allee_threshold_ssb);
      ctx.save();
      ctx.setLineDash([5, 4]);
      ctx.strokeStyle = "#c53929";
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(pad.left, yThr);
      ctx.lineTo(pad.left + pw, yThr);
      ctx.stroke();
      ctx.fillStyle = "#c53929";
      ctx.textAlign = "right";
      ctx.font = "600 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.fillText("Cultivation Threshold (" + Math.round(DATA.reference.allee_threshold_ssb) + " kt)", pad.left + pw - 4, yThr - 6);
      ctx.restore();
    }

    // Vertical milestone markers
    function drawVMarker(yr, label, color) {
      const x = xScale(yr);
      ctx.save();
      ctx.setLineDash([3, 3]);
      ctx.strokeStyle = color;
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      ctx.moveTo(x, pad.top);
      ctx.lineTo(x, pad.top + ph);
      ctx.stroke();
      ctx.fillStyle = color;
      ctx.font = "600 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(label, x + 5, pad.top + 11);
      ctx.restore();
    }
    drawVMarker(1992, "1992 Moratorium", "#9e2a2b");
    drawVMarker(DATA.policy_start_year, DATA.policy_start_year + " Policy Action", "#1d5fa7");

    // X-axis ticks
    ctx.fillStyle = "#5c7089";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (let yr = 1970; yr <= xMax; yr += 10) {
      const x = xScale(yr);
      ctx.fillText(String(yr), x, pad.top + ph + 8);
    }

    // Series helper
    function plotLine(series, color, widthPx) {
      ctx.strokeStyle = color;
      ctx.lineWidth = widthPx;
      ctx.lineJoin = "round";
      ctx.beginPath();
      for (let i = 0; i < years.length; i++) {
        const x = xScale(years[i]);
        const y = yScale(series[i]);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    plotLine(sculpinSeries, "#7c3aed", 2.0);
    plotLine(forageSeries, "#d97706", 2.4);
    plotLine(codSeries, "#1d5fa7", 3.1);
  }

  // Draw Stacked / Multi-line Yield Chart
  function drawYieldChart(scen) {
    const canvas = document.getElementById("canvas-yield");
    if (!canvas) return;
    const { ctx, width, height } = setupCanvas(canvas);
    ctx.clearRect(0, 0, width, height);

    const hist = DATA.history;
    const projYears = DATA.proj_years;
    const years = hist.years.concat(projYears.slice(1));
    const codY = hist.cod_yield.concat(scen.cod_yield.slice(1));
    const pelY = hist.pelagic_yield.concat(scen.pelagic_yield.slice(1));
    const othY = hist.other_yield.concat(scen.other_yield.slice(1));

    const pad = { left: 56, right: 20, top: 14, bottom: 28 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const xMin = years[0];
    const xMax = years[years.length - 1];
    const maxY = Math.max(250, Math.ceil(Math.max(...codY, ...pelY, ...othY) * 1.15 / 50) * 50);

    const xScale = (yr) => pad.left + ((yr - xMin) / (xMax - xMin)) * pw;
    const yScale = (v) => pad.top + ph - (v / maxY) * ph;

    ctx.strokeStyle = "#e5ecf3";
    ctx.lineWidth = 1;
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#5c7089";
    for (let i = 0; i <= 4; i++) {
      const val = (maxY / 4) * i;
      const y = yScale(val);
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(pad.left + pw, y);
      ctx.stroke();
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(Math.round(val) + " kt", pad.left - 7, y);
    }

    const xPolicy = xScale(DATA.policy_start_year);
    ctx.save();
    ctx.setLineDash([3, 3]);
    ctx.strokeStyle = "#1d5fa7";
    ctx.beginPath();
    ctx.moveTo(xPolicy, pad.top);
    ctx.lineTo(xPolicy, pad.top + ph);
    ctx.stroke();
    ctx.restore();

    ctx.fillStyle = "#5c7089";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (let yr = 1970; yr <= xMax; yr += 10) {
      ctx.fillText(String(yr), xScale(yr), pad.top + ph + 6);
    }

    function plotLine(series, color, lw) {
      ctx.strokeStyle = color;
      ctx.lineWidth = lw;
      ctx.beginPath();
      for (let i = 0; i < years.length; i++) {
        const x = xScale(years[i]);
        const y = yScale(series[i]);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
    }

    plotLine(othY, "#15805c", 1.9);
    plotLine(pelY, "#d97706", 2.3);
    plotLine(codY, "#1d5fa7", 2.6);
  }

  // Update Simulator KPIs and Charts
  function updateSimulator() {
    document.getElementById("val-fp").textContent = state.fp.toFixed(2) + " /yr";
    document.getElementById("val-fm").textContent = state.fm.toFixed(2) + " /yr";
    document.getElementById("val-fcod").textContent = state.fcod.toFixed(2) + " /yr";
    document.getElementById("val-pulse").textContent = state.taper ? "On (Yr 20)" : "Off";

    const scen = getInterpolatedScenario(state.fp, state.fm, state.fcod, state.taper);
    const n = scen.cod_ssb.length - 1;
    const endCod = scen.cod_ssb[n];
    const endForage = scen.forage_bio[n];
    const endCodY = scen.cod_yield[n];
    const endTotY = endCodY + scen.pelagic_yield[n] + scen.other_yield[n];
    const pctHist = (endCod / DATA.reference.precollapse_cod_ssb) * 100;

    document.getElementById("kpi-cod-ssb").textContent = Math.round(endCod) + " kt";
    document.getElementById("kpi-cod-pct").textContent = Math.round(pctHist) + "% of pre-collapse SSB";
    document.getElementById("kpi-forage-bio").textContent = Math.round(endForage) + " kt";
    document.getElementById("kpi-total-yield").textContent = Math.round(endTotY) + " kt/yr";
    document.getElementById("kpi-cod-yield").textContent = "Cod yield: " + Math.round(endCodY) + " kt/yr";
    document.getElementById("kpi-larval-m").textContent = scen.larval_mort_end.toFixed(2) + " /yr";

    const banner = document.getElementById("regime-status-banner");
    const statusText = document.getElementById("regime-status-text");
    const statusNote = document.getElementById("regime-threshold-note");
    banner.classList.remove("trapped", "recovered", "transitional");

    if (endCod >= DATA.reference.allee_threshold_ssb * 2.2) {
      banner.classList.add("recovered");
      statusText.textContent = "Regime Outcome: Cod Rebuilds to Self-Sustaining Cultivated State";
      statusNote.textContent = "Adult cod restores top-down forage control (" + Math.round(pctHist) + "% of historical SSB)";
    } else if (endCod >= DATA.reference.allee_threshold_ssb * 0.9) {
      banner.classList.add("transitional");
      statusText.textContent = "Regime Outcome: Near Tipping Threshold (Slow / Partial Transition)";
      statusNote.textContent = "Slightly stronger forage/mesopredator harvest or lower cod F accelerates recovery";
    } else {
      banner.classList.add("trapped");
      statusText.textContent = "Regime Outcome: Trapped in Forage-Dominated Collapsed State";
      statusNote.textContent = "Early-life predation & competition suppress juvenile cod recruitment";
    }

    drawTimeSeries(scen);
    drawYieldChart(scen);
  }

  // Draw Size-Spectrum / Mortality / Feeding Diagnostic View
  function drawSpectraDiagnostics() {
    const canvas = document.getElementById("canvas-spectra");
    if (!canvas || !DATA.size_diagnostics) return;
    const { ctx, width, height } = setupCanvas(canvas);
    ctx.clearRect(0, 0, width, height);

    const diag = DATA.size_diagnostics;
    const logW = diag.log10_w; // e.g. -3 to 4.3
    const pad = { left: 66, right: 24, top: 24, bottom: 40 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const xMin = logW[0];
    const xMax = logW[logW.length - 1];
    const xScale = (lw) => pad.left + ((lw - xMin) / (xMax - xMin)) * pw;

    // Highlight vulnerable larval/juvenile Cod nursery window (0.002g to 15g)
    const xNursery0 = xScale(Math.log10(0.002));
    const xNursery1 = xScale(Math.log10(15));
    ctx.fillStyle = "rgba(217, 119, 6, 0.08)";
    ctx.fillRect(xNursery0, pad.top, xNursery1 - xNursery0, ph);
    ctx.fillStyle = "#b45309";
    ctx.font = "600 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Vulnerable Cod Egg / Larval / 0-Group Window", (xNursery0 + xNursery1) / 2, pad.top + 12);

    const captionEl = document.getElementById("spectra-caption");

    if (state.spectrumView === "spectra") {
      // Log10 biomass per log weight bin
      const yMin = 8.2;
      const yMax = 12.6;
      const yScale = (v) => pad.top + ph - ((v - yMin) / (yMax - yMin)) * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, yMin, yMax, yScale,
        (v) => "10^" + (v - 9).toFixed(0) + " t", "Biomass per log-size octave");

      const colors = {
        Forage: "#d97706",
        Sculpin: "#7c3aed",
        Cod: "#1d5fa7"
      };

      function drawRegimeSpectra(regKey, dashed) {
        const reg = diag[regKey];
        ["Forage", "Sculpin", "Cod"].forEach((g) => {
          const arr = reg.spectra_log10[g];
          ctx.save();
          if (dashed) ctx.setLineDash([6, 4]);
          ctx.strokeStyle = colors[g];
          ctx.lineWidth = g === "Cod" ? 2.8 : 2.2;
          ctx.beginPath();
          let started = false;
          for (let i = 0; i < logW.length; i++) {
            const val = arr[i];
            if (val === null || val < yMin) {
              started = false;
              continue;
            }
            const x = xScale(logW[i]);
            const y = yScale(val);
            if (!started) {
              ctx.moveTo(x, y);
              started = true;
            } else {
              ctx.lineTo(x, y);
            }
          }
          ctx.stroke();
          ctx.restore();
        });
      }

      if (state.regimeCompare === "moratorium") {
        drawRegimeSpectra("moratorium", false);
      } else if (state.regimeCompare === "cultivated") {
        drawRegimeSpectra("cultivated", false);
      } else {
        drawRegimeSpectra("cultivated", true);
        drawRegimeSpectra("moratorium", false);
      }

      drawInlineLegend(ctx, pad.left + pw - 330, pad.top + 10, [
        { label: "Atlantic Cod", color: "#1d5fa7" },
        { label: "Pelagic Forage Guild", color: "#d97706" },
        { label: "Sculpin Mesopredator", color: "#7c3aed" }
      ], state.regimeCompare === "both" ? "(Solid = Moratorium, Dashed = Cultivated)" : "");

      captionEl.innerHTML =
        "<strong>Biomass size spectra across body weight:</strong> In the 1992 moratorium regime, the pelagic forage dome (10–600 g, orange) and sculpin dome (50–900 g, purple) expand dramatically, while the adult cod dome (1–20 kg, blue) is truncated by more than an order of magnitude.";
    } else if (state.spectrumView === "mortality") {
      const mMor = diag.moratorium.cod_pred_mort;
      const mCul = diag.cultivated.cod_pred_mort;
      const yMax = Math.ceil(Math.max(...mMor, ...mCul) * 1.15 * 2) / 2;
      const yScale = (v) => pad.top + ph - (v / yMax) * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, 0, yMax, yScale,
        (v) => v.toFixed(1) + " /yr", "Predation mortality on Cod (1/yr)");

      function plotCurve(arr, color, dashed, lw) {
        ctx.save();
        if (dashed) ctx.setLineDash([6, 4]);
        ctx.strokeStyle = color;
        ctx.lineWidth = lw;
        ctx.beginPath();
        for (let i = 0; i < logW.length; i++) {
          const x = xScale(logW[i]);
          const y = yScale(arr[i]);
          if (i === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.restore();
      }

      if (state.regimeCompare === "moratorium") {
        plotCurve(mMor, "#c53929", false, 3.0);
      } else if (state.regimeCompare === "cultivated") {
        plotCurve(mCul, "#15805c", false, 3.0);
      } else {
        plotCurve(mCul, "#15805c", true, 2.6);
        plotCurve(mMor, "#c53929", false, 3.0);
      }

      drawInlineLegend(ctx, pad.left + pw - 310, pad.top + 10, [
        { label: "Moratorium Trap (Forage-Locked)", color: "#c53929" },
        { label: "Cultivated / Recovered Regime", color: "#15805c" }
      ], "");

      captionEl.innerHTML =
        "<strong>Early-life predator pit:</strong> Abundant herring, mackerel, and sculpin elevate instantaneous predation mortality on 1 mg–20 g cod several-fold. Because cohort survival integrates exp(&minus;&int; (&mu;/<em>g</em>) d<em>w</em>), this bottleneck slashes survival from egg to maturity below replacement.";
    } else {
      const fMor = diag.moratorium.cod_feeding;
      const fCul = diag.cultivated.cod_feeding;
      const yScale = (v) => pad.top + ph - v * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, 0, 1.0, yScale,
        (v) => v.toFixed(2), "Cod feeding level f(w)");

      function plotCurve(arr, color, dashed, lw) {
        ctx.save();
        if (dashed) ctx.setLineDash([6, 4]);
        ctx.strokeStyle = color;
        ctx.lineWidth = lw;
        ctx.beginPath();
        for (let i = 0; i < logW.length; i++) {
          const x = xScale(logW[i]);
          const y = yScale(arr[i]);
          if (i === 0) ctx.moveTo(x, y);
          else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.restore();
      }

      if (state.regimeCompare === "moratorium") {
        plotCurve(fMor, "#c53929", false, 2.8);
      } else if (state.regimeCompare === "cultivated") {
        plotCurve(fCul, "#15805c", false, 2.8);
      } else {
        plotCurve(fCul, "#15805c", true, 2.5);
        plotCurve(fMor, "#c53929", false, 2.8);
      }

      drawInlineLegend(ctx, pad.left + pw - 310, pad.top + 10, [
        { label: "Moratorium Trap", color: "#c53929" },
        { label: "Cultivated Regime", color: "#15805c" }
      ], "");

      captionEl.innerHTML =
        "<strong>Dual squeeze on juvenile growth:</strong> Overabundant planktivores graze down the shared zooplankton spectrum (0.001–10 g), lowering larval/juvenile cod feeding levels and slowing growth through the high-mortality size window.";
    }
  }

  function drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, yMin, yMax, yScale, yFmt, yTitle) {
    ctx.strokeStyle = "#e5ecf3";
    ctx.lineWidth = 1;
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#5c7089";

    const ySteps = 4;
    for (let i = 0; i <= ySteps; i++) {
      const v = yMin + ((yMax - yMin) / ySteps) * i;
      const y = yScale(v);
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(pad.left + pw, y);
      ctx.stroke();
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(yFmt(v), pad.left - 8, y);
    }

    const weightTicks = [
      { lw: -3, label: "1 mg" },
      { lw: -2, label: "10 mg" },
      { lw: -1, label: "0.1 g" },
      { lw: 0, label: "1 g" },
      { lw: 1, label: "10 g" },
      { lw: 2, label: "100 g" },
      { lw: 3, label: "1 kg" },
      { lw: 4, label: "10 kg" }
    ];
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    weightTicks.forEach((t) => {
      if (t.lw >= xMin && t.lw <= xMax) {
        const x = xScale(t.lw);
        ctx.fillText(t.label, x, pad.top + ph + 8);
      }
    });

    ctx.fillStyle = "#334e68";
    ctx.font = "600 11.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillText("Individual Body Size (weight)", pad.left + pw / 2, pad.top + ph + 24);
  }

  function drawInlineLegend(ctx, x, y, items, subtitle) {
    ctx.save();
    ctx.font = "600 11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "left";
    ctx.textBaseline = "middle";
    let curY = y + 8;
    items.forEach((it) => {
      ctx.strokeStyle = it.color;
      ctx.lineWidth = 3;
      ctx.beginPath();
      ctx.moveTo(x, curY);
      ctx.lineTo(x + 18, curY);
      ctx.stroke();
      ctx.fillStyle = "#10243e";
      ctx.fillText(it.label, x + 24, curY);
      curY += 17;
    });
    if (subtitle) {
      ctx.fillStyle = "#5c7089";
      ctx.font = "10.5px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.fillText(subtitle, x, curY + 2);
    }
    ctx.restore();
  }

  // Draw Hysteresis / Bifurcation Diagrams
  function drawBifurcationCharts() {
    if (!DATA.bifurcation) return;

    // Panel A: vs F_cod
    const c1 = document.getElementById("canvas-hysteresis-cod");
    if (c1) {
      const { ctx, width, height } = setupCanvas(c1);
      ctx.clearRect(0, 0, width, height);
      const b = DATA.bifurcation.vs_fcod;
      drawBifurcationPanel(ctx, width, height, {
        xValsUp: b.f_up,
        yValsUp: b.ssb_up,
        xValsDown: b.f_down,
        yValsDown: b.ssb_down,
        xLabel: "Directed Cod Fishing Mortality F_cod (1/yr)",
        upLabel: "Upper Cultivated Branch",
        downLabel: "Collapsed Moratorium Branch",
        tippingX: b.collapse_f,
        tippingLabel: "Collapse Cliff (F ≈ " + b.collapse_f.toFixed(2) + ")",
        isEscape: false
      });
    }

    // Panel B: vs F_pelagic at F_cod = 0
    const c2 = document.getElementById("canvas-hysteresis-forage");
    if (c2) {
      const { ctx, width, height } = setupCanvas(c2);
      ctx.clearRect(0, 0, width, height);
      const b = DATA.bifurcation.vs_fpelagic;
      drawBifurcationPanel(ctx, width, height, {
        xValsUp: b.f_up,
        yValsUp: b.ssb_up,
        xValsDown: b.f_down,
        yValsDown: b.ssb_down,
        xLabel: "Pelagic Forage Harvest F_pelagic (1/yr, with F_cod = 0)",
        upLabel: "Recovered Branch (Holds on Taper)",
        downLabel: "Starting Collapsed Branch",
        tippingX: b.recovery_f,
        tippingLabel: "Recovery Threshold (F ≈ " + b.recovery_f.toFixed(2) + ")",
        isEscape: true
      });
    }
  }

  function drawBifurcationPanel(ctx, width, height, cfg) {
    const pad = { left: 54, right: 18, top: 20, bottom: 38 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const allX = cfg.xValsUp.concat(cfg.xValsDown);
    const allY = cfg.yValsUp.concat(cfg.yValsDown);
    const xMin = Math.min(...allX);
    const xMax = Math.max(...allX);
    const yMax = Math.ceil((Math.max(...allY) * 1.15) / 100) * 100;

    const xScale = (x) => pad.left + ((x - xMin) / (xMax - xMin)) * pw;
    const yScale = (y) => pad.top + ph - (y / yMax) * ph;

    ctx.strokeStyle = "#e5ecf3";
    ctx.lineWidth = 1;
    ctx.font = "11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillStyle = "#5c7089";

    for (let i = 0; i <= 4; i++) {
      const val = (yMax / 4) * i;
      const y = yScale(val);
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(pad.left + pw, y);
      ctx.stroke();
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      ctx.fillText(Math.round(val) + " kt", pad.left - 6, y);
    }

    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    const xTicks = 5;
    for (let i = 0; i <= xTicks; i++) {
      const val = xMin + ((xMax - xMin) / xTicks) * i;
      ctx.fillText(val.toFixed(2), xScale(val), pad.top + ph + 7);
    }

    ctx.fillStyle = "#334e68";
    ctx.font = "600 11px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.fillText(cfg.xLabel, pad.left + pw / 2, pad.top + ph + 23);

    // Tipping vertical line
    if (cfg.tippingX !== undefined) {
      const xt = xScale(cfg.tippingX);
      ctx.save();
      ctx.setLineDash([4, 4]);
      ctx.strokeStyle = cfg.isEscape ? "#15805c" : "#c53929";
      ctx.lineWidth = 1.4;
      ctx.beginPath();
      ctx.moveTo(xt, pad.top);
      ctx.lineTo(xt, pad.top + ph);
      ctx.stroke();
      ctx.fillStyle = cfg.isEscape ? "#15805c" : "#c53929";
      ctx.textAlign = cfg.tippingX > (xMin + xMax) * 0.55 ? "right" : "left";
      ctx.fillText(cfg.tippingLabel, xt + (cfg.tippingX > (xMin + xMax) * 0.55 ? -6 : 6), pad.top + 10);
      ctx.restore();
    }

    function drawBranch(xs, ys, color, lw, dashed) {
      ctx.save();
      if (dashed) ctx.setLineDash([6, 4]);
      ctx.strokeStyle = color;
      ctx.lineWidth = lw;
      ctx.beginPath();
      for (let i = 0; i < xs.length; i++) {
        const x = xScale(xs[i]);
        const y = yScale(ys[i]);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
      }
      ctx.stroke();
      ctx.restore();
    }

    if (!cfg.isEscape) {
      drawBranch(cfg.xValsUp, cfg.yValsUp, "#1d5fa7", 2.8, false);
      drawBranch(cfg.xValsDown, cfg.yValsDown, "#c53929", 2.6, true);
      drawInlineLegend(ctx, pad.left + 10, pad.top + ph - 48, [
        { label: cfg.upLabel, color: "#1d5fa7" },
        { label: cfg.downLabel, color: "#c53929" }
      ], "");
    } else {
      drawBranch(cfg.xValsDown, cfg.yValsDown, "#15805c", 2.6, true);
      drawBranch(cfg.xValsUp, cfg.yValsUp, "#1d5fa7", 2.8, false);
      drawInlineLegend(ctx, pad.left + 10, pad.top + 28, [
        { label: "Forward Scan from Collapsed State", color: "#1d5fa7" },
        { label: cfg.upLabel, color: "#15805c" }
      ], "");
    }
  }

  // Populate Policy Comparison Table
  function populatePolicyTable() {
    const tbody = document.getElementById("policy-comparison-tbody");
    if (!tbody || !DATA.policy_table) return;
    tbody.innerHTML = DATA.policy_table
      .map((row) => {
        const badgeCls = row.recovered ? "pass" : "fail";
        return `<tr>
          <td><strong>${row.name}</strong></td>
          <td>${row.fcod}</td>
          <td>${row.fpel}</td>
          <td>${row.fmeso}</td>
          <td><strong>${row.cod_ssb_kt} kt</strong> (${row.cod_pct}%)</td>
          <td>${row.total_yield_kt} kt/yr</td>
          <td><span class="badge-outcome ${badgeCls}">${row.outcome}</span></td>
        </tr>`;
      })
      .join("");
  }

  // Bind UI Controls
  function bindEvents() {
    const sFp = document.getElementById("slider-fp");
    const sFm = document.getElementById("slider-fm");
    const sFc = document.getElementById("slider-fcod");
    const cTaper = document.getElementById("check-taper");
    const presetBtns = document.querySelectorAll("#preset-buttons .preset-btn");

    function clearActivePreset() {
      presetBtns.forEach((b) => b.classList.remove("active"));
    }

    sFp.addEventListener("input", () => {
      state.fp = parseFloat(sFp.value);
      clearActivePreset();
      updateSimulator();
    });
    sFm.addEventListener("input", () => {
      state.fm = parseFloat(sFm.value);
      clearActivePreset();
      updateSimulator();
    });
    sFc.addEventListener("input", () => {
      state.fcod = parseFloat(sFc.value);
      clearActivePreset();
      updateSimulator();
    });
    cTaper.addEventListener("change", () => {
      state.taper = cTaper.checked;
      clearActivePreset();
      updateSimulator();
    });

    presetBtns.forEach((btn) => {
      btn.addEventListener("click", () => {
        const key = btn.getAttribute("data-preset");
        const p = PRESETS[key];
        if (!p) return;
        presetBtns.forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        state.fp = p.fp;
        state.fm = p.fm;
        state.fcod = p.fcod;
        state.taper = p.taper;
        sFp.value = p.fp;
        sFm.value = p.fm;
        sFc.value = p.fcod;
        cTaper.checked = p.taper;
        updateSimulator();
      });
    });

    document.querySelectorAll("[data-spectrum-view]").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.querySelectorAll("[data-spectrum-view]").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        state.spectrumView = btn.getAttribute("data-spectrum-view");
        drawSpectraDiagnostics();
      });
    });

    document.querySelectorAll("[data-regime-compare]").forEach((btn) => {
      btn.addEventListener("click", () => {
        document.querySelectorAll("[data-regime-compare]").forEach((b) => b.classList.remove("active"));
        btn.classList.add("active");
        state.regimeCompare = btn.getAttribute("data-regime-compare");
        drawSpectraDiagnostics();
      });
    });

    window.addEventListener("resize", () => {
      updateSimulator();
      drawSpectraDiagnostics();
      drawBifurcationCharts();
    });
  }

  // Init
  populatePolicyTable();
  bindEvents();
  updateSimulator();
  drawSpectraDiagnostics();
  drawBifurcationCharts();
})();
