(function () {
  "use strict";

  const DATA = window.MIZER_MODEL_DATA;
  if (!DATA) {
    console.error("MIZER_MODEL_DATA not found. Run R/03_generate_figures_and_web_data.R first.");
    return;
  }

  const FP_AX   = DATA.grid_axes.fp;   // [0.00, 0.15, 0.35, 0.50, 0.55, 0.65, 0.80]
  const FM_AX   = DATA.grid_axes.fm;   // [0.00, 0.05, 0.20, 0.35, 0.45, 0.60]
  const FCOD_AX = DATA.grid_axes.fcod; // [0.00, 0.10, 0.15, 0.25, 0.40]

  function nearestIdx(arr, val) {
    let best = 0;
    let bestDiff = Infinity;
    for (let i = 0; i < arr.length; i++) {
      const d = Math.abs(arr[i] - val);
      if (d < bestDiff) {
        bestDiff = d;
        best = i;
      }
    }
    return best;
  }

  // State stores exact grid indices — zero cross-basin interpolation
  const state = {
    fpIdx: nearestIdx(FP_AX, 0.15),
    fmIdx: nearestIdx(FM_AX, 0.05),
    fcIdx: nearestIdx(FCOD_AX, 0.00),
    taper: false,
    activePreset: "status_quo",
    spectrumView: "spectra",
    regimeCompare: "moratorium"
  };

  const PRESETS = DATA.presets || {
    status_quo:        { fp: 0.15, fm: 0.05, fcod: 0.00, taper: false },
    pelagic_rebalance: { fp: 0.55, fm: 0.05, fcod: 0.00, taper: false },
    ecosystem_pulse:   { fp: 0.65, fm: 0.45, fcod: 0.15, taper: true  },
    premature_reopen:  { fp: 0.15, fm: 0.05, fcod: 0.25, taper: false }
  };

  if (DATA.hero) {
    const setTxt = (id, val) => {
      const el = document.getElementById(id);
      if (el && val !== undefined) el.innerHTML = val;
    };
    setTxt("hero-high-low", DATA.hero.high_low_ssb);
    setTxt("hero-forage-surge", DATA.hero.forage_surge);
    setTxt("hero-survival-drop", DATA.hero.survival_drop);
    setTxt("hero-bycatch-limit", DATA.hero.bycatch_limit);
  }

  // Exact grid node lookup (no interpolation)
  function getExactScenario(i, j, k, taper) {
    const grid = taper ? DATA.grid_taper : DATA.grid_sustained;
    return grid[`${i}_${j}_${k}`];
  }

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

  // Draw Time Series Chart in Relative Years (-61..0 collapse+closure, 0..+40 intervention)
  function drawTimeSeries(scen) {
    const canvas = document.getElementById("canvas-timeseries");
    if (!canvas) return;
    const { ctx, width, height } = setupCanvas(canvas);
    ctx.clearRect(0, 0, width, height);

    const hist = DATA.history;
    const projYears = DATA.proj_years; // 0..40
    const years = hist.years.concat(projYears.slice(1));
    const codSeries = hist.cod_ssb.concat(scen.cod_ssb.slice(1));
    const forageSeries = hist.forage_bio.concat(scen.forage_bio.slice(1));
    const sculpinSeries = hist.sculpin_bio.concat(scen.sculpin_bio.slice(1));

    const pad = { left: 58, right: 20, top: 22, bottom: 36 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const xMin = years[0]; // -61
    const xMax = years[years.length - 1]; // +40
    const rawMax = Math.max(...codSeries, ...forageSeries, DATA.reference.baseline_f15_cod_ssb * 1.15);
    const yMax = Math.ceil(rawMax / 200) * 200;

    const xScale = (yr) => pad.left + ((yr - xMin) / (xMax - xMin)) * pw;
    const yScale = (v) => pad.top + ph - (v / yMax) * ph;

    const xClosure = xScale(DATA.closure_year);       // -34
    const xPolicy  = xScale(DATA.policy_start_year);  // 0
    const xSwitch  = xScale(20);                      // +20

    // Closure window shading (-34 -> 0)
    ctx.fillStyle = "rgba(197, 57, 41, 0.055)";
    ctx.fillRect(xClosure, pad.top, xPolicy - xClosure, ph);

    // Intervention window shading (0 -> +40)
    ctx.fillStyle = "rgba(29, 95, 167, 0.05)";
    ctx.fillRect(xPolicy, pad.top, pad.left + pw - xPolicy, ph);

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
      ctx.fillText(label, x + 4, pad.top + 11);
      ctx.restore();
    }
    drawVMarker(DATA.closure_year, "Closure (-34 yr)", "#9e2a2b");
    drawVMarker(DATA.policy_start_year, "Intervention (0 yr)", "#1d5fa7");
    drawVMarker(20, "Switch (+20 yr)", "#475569");

    ctx.fillStyle = "#5c7089";
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (let yr = -60; yr <= xMax; yr += 20) {
      const x = xScale(yr);
      const lbl = yr > 0 ? "+" + yr : String(yr);
      ctx.fillText(lbl + " yr", x, pad.top + ph + 8);
    }

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

    const pad = { left: 58, right: 20, top: 14, bottom: 28 };
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
    for (let yr = -60; yr <= xMax; yr += 20) {
      const lbl = yr > 0 ? "+" + yr : String(yr);
      ctx.fillText(lbl + " yr", xScale(yr), pad.top + ph + 6);
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

  function updateSimulator() {
    const fpVal   = FP_AX[state.fpIdx];
    const fmVal   = FM_AX[state.fmIdx];
    const fcodVal = FCOD_AX[state.fcIdx];

    document.getElementById("val-fp").textContent   = fpVal.toFixed(2) + " /yr";
    document.getElementById("val-fm").textContent   = fmVal.toFixed(2) + " /yr";
    document.getElementById("val-fcod").textContent = fcodVal.toFixed(2) + " /yr";
    document.getElementById("val-pulse").textContent = state.taper ? "On (tau = 20)" : "Off";

    const scen = getExactScenario(state.fpIdx, state.fmIdx, state.fcIdx, state.taper);
    const n = scen.cod_ssb.length - 1;
    const endCod = scen.t40_cod_ssb;
    const endForage = scen.forage_bio[n];
    const endCodY = scen.cod_yield[n];
    const endTotY = endCodY + scen.pelagic_yield[n] + scen.other_yield[n];

    document.getElementById("kpi-cod-ssb").textContent = endCod.toFixed(1) + " kt";
    document.getElementById("kpi-cod-asymp").textContent = "Asymptotic (tau=340): " + scen.asymp_cod_ssb.toFixed(1) + " kt";
    document.getElementById("kpi-forage-bio").textContent = Math.round(endForage) + " kt";
    document.getElementById("kpi-min-capelin").textContent = "Min Capelin (0-40 yr): " + Math.round(scen.min_capelin_bio) + " kt";
    document.getElementById("kpi-total-yield").textContent = Math.round(endTotY) + " kt/yr";
    document.getElementById("kpi-cum-yield").textContent =
      "40-yr cum.: " + Math.round(scen.cum40_tot_yield) + " kt (Cod " + Math.round(scen.cum40_cod_yield) + " kt)";
    document.getElementById("kpi-larval-m").textContent = scen.larval_mort_end.toFixed(2) + " /yr";

    // Drive banner strictly from directly simulated asymptotic outcome (NO scalar SSB multipliers)
    const banner = document.getElementById("regime-status-banner");
    const statusText = document.getElementById("regime-status-text");
    const statusNote = document.getElementById("regime-threshold-note");
    banner.classList.remove("trapped", "recovered", "transitional");

    if (scen.asymp_recovered && endCod >= 300) {
      banner.classList.add("recovered");
      statusText.textContent = "Verified Asymptotic Outcome: Upper Cultivated Basin (" + scen.asymp_cod_ssb.toFixed(0) + " kt)";
      statusNote.textContent = "Exact grid simulation: tau = 20 SSB = " + scen.t20_cod_ssb.toFixed(1) +
        " kt; tau = 40 SSB = " + endCod.toFixed(1) + " kt; settles to " + scen.asymp_cod_ssb.toFixed(1) + " kt";
    } else if (scen.asymp_recovered && endCod < 300) {
      banner.classList.add("transitional");
      statusText.textContent = "Verified Asymptotic Outcome: Slow Multi-Century Escape (" + scen.asymp_cod_ssb.toFixed(0) + " kt)";
      statusNote.textContent = "Only " + endCod.toFixed(1) + " kt by tau = 40 yr, but full-state continuation reaches " +
        scen.asymp_cod_ssb.toFixed(1) + " kt by tau = 340 yr";
    } else {
      banner.classList.add("trapped");
      statusText.textContent = "Verified Asymptotic Outcome: Lower Depensatory Basin (" + scen.asymp_cod_ssb.toFixed(1) + " kt)";
      statusNote.textContent = "Post-switch ecosystem state settles to " + scen.asymp_cod_ssb.toFixed(1) +
        " kt Cod SSB under constant final fishing efforts";
    }

    drawTimeSeries(scen);
    drawYieldChart(scen);
  }

  // Static Reference Comparison (Section 3)
  function drawSpectraDiagnostics() {
    const canvas = document.getElementById("canvas-spectra");
    if (!canvas || !DATA.size_diagnostics) return;
    const { ctx, width, height } = setupCanvas(canvas);
    ctx.clearRect(0, 0, width, height);

    const diag = DATA.size_diagnostics;
    const logW = diag.log10_w;
    const pad = { left: 66, right: 24, top: 24, bottom: 40 };
    const pw = width - pad.left - pad.right;
    const ph = height - pad.top - pad.bottom;

    const xMin = logW[0];
    const xMax = logW[logW.length - 1];
    const xScale = (lw) => pad.left + ((lw - xMin) / (xMax - xMin)) * pw;

    const xNursery0 = xScale(Math.log10(0.001));
    const xNursery1 = xScale(Math.log10(15));
    ctx.fillStyle = "rgba(217, 119, 6, 0.08)";
    ctx.fillRect(xNursery0, pad.top, xNursery1 - xNursery0, ph);
    ctx.fillStyle = "#b45309";
    ctx.font = "600 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Early-Life Cod Bottleneck Window (1 mg – 15 g)", (xNursery0 + xNursery1) / 2, pad.top + 12);

    const captionEl = document.getElementById("spectra-caption");

    if (state.spectrumView === "spectra") {
      const yMin = 8.2;
      const yMax = 12.6;
      const yScale = (v) => pad.top + ph - ((v - yMin) / (yMax - yMin)) * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, yMin, yMax, yScale,
        (v) => "10^" + (v - 9).toFixed(0) + " t");

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

      drawInlineLegend(ctx, pad.left + pw - 350, pad.top + 10, [
        { label: "Cod", color: "#1d5fa7" },
        { label: "Pelagic Forage Guild", color: "#d97706" },
        { label: "Sculpin", color: "#7c3aed" }
      ], state.regimeCompare === "both" ? "(Solid = Low-Cod, Dashed = High-Cod; static E_cod = 0)" : "(Static E_cod = 0 equilibrium)");

      captionEl.innerHTML =
        "<strong>Static E<sub>cod</sub> = 0 equilibrium spectra:</strong> Comparing the two coexisting zero-cod-fishing steady states (baseline non-cod fishing active), forage (10–600 g) and sculpin (50–900 g) domes expand strongly in the depensatory state while adult cod (> 1 kg) drops from 1018.2 kt to 15.8 kt SSB.";
    } else if (state.spectrumView === "mortality") {
      const mMor = diag.moratorium.cod_pred_mort;
      const mCul = diag.cultivated.cod_pred_mort;
      const yMax = Math.ceil(Math.max(...mMor, ...mCul) * 1.15 * 2) / 2;
      const yScale = (v) => pad.top + ph - (v / yMax) * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, 0, yMax, yScale,
        (v) => v.toFixed(1) + " /yr");

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

      drawInlineLegend(ctx, pad.left + pw - 330, pad.top + 10, [
        { label: "Depensatory Low-Cod (E_cod = 0)", color: "#c53929" },
        { label: "Cultivated High-Cod (E_cod = 0)", color: "#15805c" }
      ], "(Static reference comparison)");

      captionEl.innerHTML =
        "<strong>Early-life predation bottleneck (static E<sub>cod</sub> = 0 comparison):</strong> Mackerel, herring, and sculpin account for > 95% of predation hazard on cod below 15 g (capelin < 5%), raising cumulative hazard H(15 g) = &int;(&mu;/<em>g</em>) d<em>w</em> from 0.19 to 2.82 and cutting stationary cohort survival 14-fold.";
    } else {
      const fMor = diag.moratorium.cod_feeding;
      const fCul = diag.cultivated.cod_feeding;
      const yScale = (v) => pad.top + ph - v * ph;

      drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, 0, 1.0, yScale,
        (v) => v.toFixed(2));

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

      drawInlineLegend(ctx, pad.left + pw - 330, pad.top + 10, [
        { label: "Depensatory Low-Cod (E_cod = 0)", color: "#c53929" },
        { label: "Cultivated High-Cod (E_cod = 0)", color: "#15805c" }
      ], "(Static reference comparison)");

      captionEl.innerHTML =
        "<strong>Opposing ontogenetic feeding shifts (static E<sub>cod</sub> = 0 comparison):</strong> In the depensatory state, cod feeding level <em>f</em>(<em>w</em>) is <strong>1.4%–11.7% lower</strong> over 0.05–5 g (up to 15.0% lower over 0.002–15 g) due to shared resource grazing, yet <strong>higher above 200 g</strong> where surviving adult cod feed on abundant forage fish.";
    }
  }

  function drawAxesLogW(ctx, pad, pw, ph, xMin, xMax, xScale, yMin, yMax, yScale, yFmt) {
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
    ctx.fillText("Individual Body Weight w", pad.left + pw / 2, pad.top + ph + 24);
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

  function drawBifurcationCharts() {
    if (!DATA.bifurcation) return;

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
        xLabel: "Directed Cod Fishing Effort E_cod (1/yr)",
        upLabel: "Forward (from cultivated state)",
        downLabel: "Backward (from depensatory state)",
        bracketLo: b.bracket_lo,
        bracketHi: b.bracket_hi,
        bracketLabel: "Transition [" + b.bracket_lo.toFixed(3) + ", " + b.bracket_hi.toFixed(3) + "]",
        isEscape: false
      });
    }

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
        xLabel: "Pelagic Effort E_pel (1/yr, at E_cod = 0)",
        upLabel: "Backward (holds on return)",
        downLabel: "Forward (from low state)",
        bracketLo: b.bracket_lo,
        bracketHi: b.bracket_hi,
        bracketLabel: "Escape [" + b.bracket_lo.toFixed(3) + ", " + b.bracket_hi.toFixed(3) + "]",
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

    // Shaded transition bracket
    if (cfg.bracketLo !== undefined && cfg.bracketHi !== undefined) {
      const x0 = xScale(cfg.bracketLo);
      const x1 = xScale(cfg.bracketHi);
      ctx.save();
      ctx.fillStyle = cfg.isEscape ? "rgba(21, 128, 92, 0.18)" : "rgba(197, 57, 41, 0.18)";
      ctx.fillRect(x0, pad.top, Math.max(3, x1 - x0), ph);
      ctx.fillStyle = cfg.isEscape ? "#15805c" : "#c53929";
      ctx.font = "600 10.5px -apple-system, BlinkMacSystemFont, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(cfg.bracketLabel, x1 + 5, pad.top + 11);
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
      drawBranch(cfg.xValsUp, cfg.yValsUp, "#d97706", 2.8, false);
      drawInlineLegend(ctx, pad.left + 10, pad.top + 28, [
        { label: cfg.downLabel, color: "#d97706" },
        { label: cfg.upLabel, color: "#15805c" }
      ], "");
    }
  }

  function populatePolicyTable() {
    const tbody = document.getElementById("policy-comparison-tbody");
    if (!tbody || !DATA.policy_table) return;
    tbody.innerHTML = DATA.policy_table
      .map((row) => {
        const badgeCls = row.recovered ? "pass" : "fail";
        return `<tr>
          <td><strong>${row.name}</strong></td>
          <td>${row.ecap_ph1} / ${row.ehm_ph1} / ${row.emeso_ph1}</td>
          <td>${row.ecod_ph2}</td>
          <td>${row.tapered}</td>
          <td><strong>${row.t40_ssb_kt} kt</strong> &rarr; ${row.asymp_ssb_kt} kt</td>
          <td>${row.cum40_tot_yld_kt} kt &middot; (${row.t40_tot_yld_kt} kt/yr)</td>
          <td>${row.min_cap_ph1_kt} kt</td>
          <td><span class="badge-outcome ${badgeCls}">${row.outcome}</span></td>
        </tr>`;
      })
      .join("");
  }

  function bindEvents() {
    const sFp = document.getElementById("slider-fp");
    const sFm = document.getElementById("slider-fm");
    const sFc = document.getElementById("slider-fcod");
    const cTaper = document.getElementById("check-taper");
    const presetBtns = document.querySelectorAll("#preset-buttons .preset-btn");

    sFp.max = String(FP_AX.length - 1);
    sFm.max = String(FM_AX.length - 1);
    sFc.max = String(FCOD_AX.length - 1);
    sFp.value = String(state.fpIdx);
    sFm.value = String(state.fmIdx);
    sFc.value = String(state.fcIdx);

    function clearActivePreset() {
      presetBtns.forEach((b) => b.classList.remove("active"));
    }

    sFp.addEventListener("input", () => {
      state.fpIdx = parseInt(sFp.value, 10);
      clearActivePreset();
      updateSimulator();
    });
    sFm.addEventListener("input", () => {
      state.fmIdx = parseInt(sFm.value, 10);
      clearActivePreset();
      updateSimulator();
    });
    sFc.addEventListener("input", () => {
      state.fcIdx = parseInt(sFc.value, 10);
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
        state.fpIdx = nearestIdx(FP_AX, p.fp);
        state.fmIdx = nearestIdx(FM_AX, p.fm);
        state.fcIdx = nearestIdx(FCOD_AX, p.fcod);
        state.taper = p.taper;
        sFp.value = String(state.fpIdx);
        sFm.value = String(state.fmIdx);
        sFc.value = String(state.fcIdx);
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

  populatePolicyTable();
  bindEvents();
  updateSimulator();
  drawSpectraDiagnostics();
  drawBifurcationCharts();
})();
