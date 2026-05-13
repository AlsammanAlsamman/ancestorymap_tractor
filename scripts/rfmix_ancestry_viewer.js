(function () {
  const state = {
    analysis: null,
    qFiles: [],
    qSummary: null,
    msp: null,
    loci: [],
    genes: [],
  };

  const defaultFiles = {
    analysis: "../configs/analysis.yml",
    qFiles: ["../results_hispanic/rfmix/chr1.deconvoluted.rfmix.Q"],
    msp: "../results_hispanic/rfmix/chr1.deconvoluted.msp.tsv",
    genes: "../resources/NCBI37.3.gene.loc",
    loci: "../input/loci.txt",
  };

  const ui = {};

  document.addEventListener("DOMContentLoaded", init);

  function init() {
    bindUi();
    setStatus("Choose files or load defaults.");
  }

  function bindUi() {
    ui.analysisFile = document.getElementById("analysisFile");
    ui.qFiles = document.getElementById("qFiles");
    ui.mspFile = document.getElementById("mspFile");
    ui.geneFile = document.getElementById("geneFile");
    ui.lociFile = document.getElementById("lociFile");
    ui.chromosomeSelect = document.getElementById("chromosomeSelect");
    ui.locusSelect = document.getElementById("locusSelect");
    ui.sortAncestrySelect = document.getElementById("sortAncestrySelect");
    ui.loadDefaultsBtn = document.getElementById("loadDefaultsBtn");
    ui.renderBtn = document.getElementById("renderBtn");
    ui.status = document.getElementById("status");
    ui.summaryGrid = document.getElementById("summaryGrid");
    ui.ancestryLegend = document.getElementById("ancestryLegend");
    ui.regionSvg = document.getElementById("regionSvg");
    ui.globalContributionPlot = document.getElementById("globalContributionPlot");
    ui.locusContributionPlot = document.getElementById("locusContributionPlot");
    ui.sampleHeatmapPlot = document.getElementById("sampleHeatmapPlot");
    ui.sampleStackedPlot = document.getElementById("sampleStackedPlot");
    ui.locusTableBody = document.getElementById("locusTableBody");

    ui.loadDefaultsBtn.addEventListener("click", loadDefaultFiles);
    ui.renderBtn.addEventListener("click", loadAndRenderFromInputs);
    ui.chromosomeSelect.addEventListener("change", handleChromosomeChange);
    ui.locusSelect.addEventListener("change", renderAll);
    ui.sortAncestrySelect.addEventListener("change", renderSamplePlots);
  }

  function handleChromosomeChange() {
    updateLocusOptions();
    renderAll();
  }

  async function loadDefaultFiles() {
    try {
      setStatus("Loading default repository files...");
      const [analysisText, qTexts, mspText, geneText, lociText] = await Promise.all([
        fetchText(defaultFiles.analysis),
        Promise.all(defaultFiles.qFiles.map((path) => fetchText(path).then((text) => ({ path, text })))),
        fetchText(defaultFiles.msp),
        fetchText(defaultFiles.genes),
        fetchText(defaultFiles.loci),
      ]);

      state.analysis = parseAnalysis(analysisText);
      state.qFiles = qTexts.map(({ path, text }) => parseQFile(text, path.split("/").pop()));
      state.qSummary = summarizeQFiles(state.qFiles, state.analysis);
      state.msp = parseMsp(mspText, state.analysis);
      state.genes = parseGeneLoc(geneText);
      state.loci = parseLoci(lociText);
      populateControls();
      renderAll();
      setStatus("Default repository files loaded.");
    } catch (error) {
      setStatus(`Failed to load default files: ${error.message}`, true);
    }
  }

  async function loadAndRenderFromInputs() {
    try {
      setStatus("Reading selected files...");
      const analysisFile = requireFile(ui.analysisFile.files[0], "analysis.yml");
      const qFileList = Array.from(ui.qFiles.files || []);
      const mspFile = requireFile(ui.mspFile.files[0], "RFMix MSP file");
      const geneFile = requireFile(ui.geneFile.files[0], "gene annotation file");
      const lociFile = requireFile(ui.lociFile.files[0], "loci file");
      if (!qFileList.length) {
        throw new Error("Please choose at least one RFMix Q file.");
      }

      const [analysisText, qTexts, mspText, geneText, lociText] = await Promise.all([
        readFileText(analysisFile),
        Promise.all(qFileList.map(async (file) => ({ name: file.name, text: await readFileText(file) }))),
        readFileText(mspFile),
        readFileText(geneFile),
        readFileText(lociFile),
      ]);

      state.analysis = parseAnalysis(analysisText);
      state.qFiles = qTexts.map(({ name, text }) => parseQFile(text, name));
      state.qSummary = summarizeQFiles(state.qFiles, state.analysis);
      state.msp = parseMsp(mspText, state.analysis);
      state.genes = parseGeneLoc(geneText);
      state.loci = parseLoci(lociText);
      populateControls();
      renderAll();
      setStatus("Viewer updated from selected files.");
    } catch (error) {
      setStatus(error.message, true);
    }
  }

  function parseAnalysis(text) {
    const yaml = window.jsyaml.load(text);
    if (!yaml || !yaml.reference || !yaml.ancestry) {
      throw new Error("analysis.yml is missing required reference/ancestry sections.");
    }

    const populations = yaml.reference.populations || [];
    const popToLabel = yaml.ancestry.population_to_label || {};
    const labelToColor = yaml.ancestry.label_to_color || {};
    if (!Array.isArray(populations) || !populations.length) {
      throw new Error("analysis.yml reference.populations must be a non-empty list.");
    }

    const ancestries = populations.map((population, index) => {
      const label = popToLabel[population];
      if (!label) {
        throw new Error(`analysis.yml missing ancestry label for population ${population}`);
      }
      const color = labelToColor[label];
      if (!color) {
        throw new Error(`analysis.yml missing ancestry color for label ${label}`);
      }
      return {
        population,
        label,
        code: index,
        color,
      };
    });

    return {
      yaml,
      ancestries,
      codeToAncestry: new Map(ancestries.map((entry) => [String(entry.code), entry])),
      labelToAncestry: new Map(ancestries.map((entry) => [entry.label, entry])),
    };
  }

  function parseQFile(text, sourceName) {
    const lines = text.split(/\r?\n/).filter(Boolean);
    const headerLine = lines.find((line) => line.startsWith("#sample") || line.startsWith("sample\t") || line.startsWith("sample "));
    if (!headerLine) {
      throw new Error(`Q file ${sourceName} is missing a #sample header.`);
    }
    const header = tokenize(headerLine.replace(/^#/, ""));
    const codeColumns = header.slice(1);
    const rows = lines
      .filter((line) => !line.startsWith("#"))
      .map((line) => tokenize(line))
      .filter((parts) => parts.length >= header.length)
      .map((parts) => ({
        sample: parts[0],
        values: codeColumns.map((code, idx) => ({ code: String(code), value: Number(parts[idx + 1] || 0) })),
      }));

    return { sourceName, codeColumns: codeColumns.map(String), rows };
  }

  function summarizeQFiles(qFiles, analysis) {
    const sampleMap = new Map();
    qFiles.forEach((qFile) => {
      qFile.rows.forEach((row) => {
        if (!sampleMap.has(row.sample)) {
          sampleMap.set(
            row.sample,
            new Map(analysis.ancestries.map((entry) => [String(entry.code), { sum: 0, count: 0 }]))
          );
        }
        const entryMap = sampleMap.get(row.sample);
        row.values.forEach(({ code, value }) => {
          if (!entryMap.has(code)) {
            entryMap.set(code, { sum: 0, count: 0 });
          }
          const slot = entryMap.get(code);
          slot.sum += value;
          slot.count += 1;
        });
      });
    });

    const samples = Array.from(sampleMap.entries()).map(([sample, codeMap]) => {
      const contributions = analysis.ancestries.map((entry) => {
        const slot = codeMap.get(String(entry.code)) || { sum: 0, count: 0 };
        return {
          ...entry,
          value: slot.count ? slot.sum / slot.count : 0,
        };
      });
      return { sample, contributions };
    });

    const global = analysis.ancestries.map((entry) => {
      const values = samples.map((sample) => sample.contributions.find((item) => item.code === entry.code)?.value || 0);
      return {
        ...entry,
        value: values.length ? mean(values) : 0,
      };
    });

    return { samples, global };
  }

  function parseMsp(text, analysis) {
    const lines = text.split(/\r?\n/).filter(Boolean);
    const headerLine = lines.find((line) => line.startsWith("#chm"));
    const dataLines = lines.filter((line) => !line.startsWith("#"));
    if (!dataLines.length) {
      throw new Error("MSP file has no ancestry rows.");
    }
    const firstParts = tokenize(dataLines[0]);
    if (firstParts.length < 8) {
      throw new Error("MSP file does not have enough columns to parse ancestry calls.");
    }
    const sampleColumns = headerLine
      ? tokenize(headerLine.replace(/^#/, "")).slice(7)
      : Array.from({ length: Math.max(0, firstParts.length - 7) }, (_, index) => `hap${index}`);
    const rows = lines
      .filter((line) => !line.startsWith("#"))
      .map((line) => tokenize(line))
      .filter((parts) => parts.length >= 8)
      .map((parts) => {
        const chrom = String(parts[0]).replace(/^chr/i, "");
        const start = Number(parts[1]);
        const end = Number(parts[2]);
        const calls = parts.slice(7);
        const counts = new Map(analysis.ancestries.map((entry) => [String(entry.code), 0]));
        calls.forEach((call) => {
          const normalized = String(call);
          if (counts.has(normalized)) {
            counts.set(normalized, counts.get(normalized) + 1);
          }
        });
        const total = calls.length || 1;
        return {
          chrom,
          start,
          end,
          windowCount: Number(parts[5]) || 0,
          snpCount: Number(parts[6]) || 0,
          calls,
          sampleColumns,
          proportions: analysis.ancestries.map((entry) => ({
            ...entry,
            value: (counts.get(String(entry.code)) || 0) / total,
          })),
        };
      });

    return {
      chromosome: rows[0].chrom,
      sampleColumns,
      rows,
    };
  }

  function parseLoci(text) {
    const lines = text.split(/\r?\n/).filter(Boolean);
    const dataLines = /^locus\s+/i.test(lines[0] || "") ? lines.slice(1) : lines;
    return dataLines
      .map((line) => tokenize(line))
      .filter((parts) => parts.length >= 4)
      .map((parts) => ({
        locus: parts[0],
        chrom: String(parts[1]).replace(/^chr/i, ""),
        start: Number(parts[2]),
        end: Number(parts[3]),
      }));
  }

  function parseGeneLoc(text) {
    return text
      .split(/\r?\n/)
      .filter(Boolean)
      .map((line) => tokenize(line))
      .filter((parts) => parts.length >= 6)
      .map((parts) => ({
        geneId: parts[0],
        chrom: String(parts[1]).replace(/^chr/i, ""),
        start: Number(parts[2]),
        end: Number(parts[3]),
        strand: parts[4],
        symbol: parts[5],
      }));
  }

  function populateControls() {
    const chromosomes = Array.from(
      new Set([
        ...state.loci.map((entry) => entry.chrom),
        ...state.genes.map((entry) => entry.chrom),
        state.msp?.chromosome,
      ].filter(Boolean))
    ).sort((a, b) => chromOrder(a) - chromOrder(b));

    fillSelect(
      ui.chromosomeSelect,
      chromosomes.map((chrom) => ({ value: chrom, label: `chr${chrom}` })),
      state.msp ? state.msp.chromosome : chromosomes[0]
    );

    updateLocusOptions();

    fillSelect(
      ui.sortAncestrySelect,
      state.analysis.ancestries.map((entry) => ({
        value: String(entry.code),
        label: `${entry.label} (code ${entry.code})`,
      })),
      String(state.analysis.ancestries[0].code)
    );
  }

  function renderAll() {
    if (!state.analysis) {
      return;
    }
    renderLegend();
    renderSummary();
    renderRegionPlot();
    renderLocusTable();
    renderContributionPlots();
    renderSamplePlots();
  }

  function renderSummary() {
    const chrom = selectedChromosome();
    const locus = selectedLocus();
    const loci = lociForChromosome(chrom);
    const genes = genesForChromosome(chrom);
    const locusGenes = locus ? genes.filter((gene) => overlaps(gene.start, gene.end, locus.start, locus.end)) : [];
    const locusGeneSymbols = unique(locusGenes.map((gene) => gene.symbol));
    const mspRows = locus ? filterRowsWithinSingleLocus(mspRowsForChromosome(chrom), locus) : [];
    const sampleCount = state.qSummary ? state.qSummary.samples.length : 0;
    const locusSize = locus ? locus.end - locus.start : 0;

    const metrics = [
      { label: "Chromosome", value: `chr${chrom}` },
      { label: "Selected locus", value: locus ? locus.locus : "-" },
      { label: "Locus size (bp)", value: locus ? formatNumber(locusSize) : "-" },
      { label: "Loaded Q files", value: String(state.qFiles.length) },
      { label: "Samples", value: String(sampleCount) },
      { label: "Loci on chromosome", value: String(loci.length) },
      { label: "Genes in selected locus", value: String(locusGeneSymbols.length) },
      { label: "MSP windows in locus", value: String(mspRows.length) },
      {
        label: "Gene names in selected locus",
        value: locusGeneSymbols.length ? locusGeneSymbols.join(", ") : "-",
        wide: true,
        long: true,
      },
    ];

    ui.summaryGrid.innerHTML = metrics
      .map(
        (metric) =>
          `<div class="metric${metric.wide ? " wide" : ""}"><div class="label">${escapeHtml(metric.label)}</div><div class="value${metric.long ? " long" : ""}">${escapeHtml(metric.value)}</div></div>`
      )
      .join("");
  }

  function renderLegend() {
    ui.ancestryLegend.innerHTML = state.analysis.ancestries
      .map(
        (entry) =>
          `<span class="legend-item"><span class="dot" style="background:${entry.color}"></span>${escapeHtml(entry.label)} (code ${entry.code})</span>`
      )
      .join("");
  }

  function renderRegionPlot() {
    const chrom = selectedChromosome();
    const loci = lociForChromosome(chrom);
    const locus = selectedLocus();
    const svg = ui.regionSvg;
    clearNode(svg);

    if (!state.msp) {
      drawSvgMessage(svg, "Load an MSP file to render locus ancestry segments.");
      return;
    }

    if (!loci.length) {
      drawSvgMessage(svg, `No loci found for chr${chrom}.`);
      return;
    }

    if (!locus) {
      drawSvgMessage(svg, `No locus selected for chr${chrom}.`);
      return;
    }

    const genes = genesForChromosome(chrom).filter((gene) => overlaps(gene.start, gene.end, locus.start, locus.end));
    const rows = filterRowsWithinSingleLocus(mspRowsForChromosome(chrom), locus);
    const locusWidth = Math.max(1, locus.end - locus.start);
    const flank = Math.max(25000, Math.round(locusWidth * 0.2));
    const minPos = Math.max(0, locus.start - flank);
    const maxPos = locus.end + flank;

    const width = 1400;
    const leftPad = 70;
    const rightPad = 20;
    const usableWidth = width - leftPad - rightPad;
    const overviewTop = 66;
    const overviewHeight = 28;
    const mspTop = 150;
    const mspHeight = 90;
    const geneTop = 300;
    const geneLanes = layoutGenes(genes);
    const laneHeight = 28;
    const totalHeight = Math.max(420, geneTop + geneLanes.laneCount * laneHeight + 60);
    svg.setAttribute("viewBox", `0 0 ${width} ${totalHeight}`);

    addSvgText(svg, leftPad, 26, `chr${chrom} locus ${locus.locus}`, "18", "700", "#18212f");
    addSvgText(svg, leftPad, 46, `${formatNumber(locus.start)}-${formatNumber(locus.end)} | ${genes.length} genes | ${rows.length} overlapping MSP windows`, "12", "400", "#5f6b7a");

    loci.forEach((locus, index) => {
      const overviewMin = Math.min(...loci.map((item) => item.start));
      const overviewMax = Math.max(...loci.map((item) => item.end));
      const x = scalePosition(locus.start, overviewMin, overviewMax, leftPad, usableWidth);
      const x2 = scalePosition(locus.end, overviewMin, overviewMax, leftPad, usableWidth);
      const selected = locus.locus === selectedLocus().locus;
      addSvgRect(
        svg,
        x,
        overviewTop,
        Math.max(3, x2 - x),
        overviewHeight,
        selected ? "rgba(15,94,168,0.32)" : index % 2 === 0 ? "rgba(15,94,168,0.10)" : "rgba(29,125,167,0.14)",
        selected ? "#0f5ea8" : "rgba(15,94,168,0.18)"
      );
      if (selected) {
        addSvgText(svg, x, overviewTop - 6, locus.locus, "11", "700", "#0f5ea8");
      }
    });

    addSvgText(svg, leftPad, overviewTop - 18, "Chromosome locus selector overview", "13", "700", "#18212f");

    rows.forEach((row) => {
      const segments = clipRowToLoci(row, [locus]);
      segments.forEach((segment) => {
        const x1 = scalePosition(segment.start, minPos, maxPos, leftPad, usableWidth);
        const x2 = scalePosition(segment.end, minPos, maxPos, leftPad, usableWidth);
        const widthPx = Math.max(1.5, x2 - x1);
        let yOffset = 0;
        state.analysis.ancestries.forEach((entry) => {
          const prop = row.proportions.find((item) => item.code === entry.code)?.value || 0;
          const height = mspHeight * prop;
          if (height > 0) {
            addSvgRect(svg, x1, mspTop + mspHeight - yOffset - height, widthPx, Math.max(1, height), entry.color, "none", 0.92);
            yOffset += height;
          }
        });
      });
    });

    const locusX1 = scalePosition(locus.start, minPos, maxPos, leftPad, usableWidth);
    const locusX2 = scalePosition(locus.end, minPos, maxPos, leftPad, usableWidth);
    addSvgRect(svg, locusX1, mspTop - 18, Math.max(2, locusX2 - locusX1), mspHeight + 36, "rgba(15,94,168,0.05)", "#0f5ea8", 1);

    addSvgLine(svg, leftPad, mspTop + mspHeight + 4, width - rightPad, mspTop + mspHeight + 4, "#94a3b8", 1);
    drawAxis(svg, minPos, maxPos, leftPad, usableWidth, mspTop + mspHeight + 24);

    geneLanes.lanes.forEach((laneGenes, laneIndex) => {
      const y = geneTop + laneIndex * laneHeight;
      laneGenes.forEach((gene) => {
        const x1 = scalePosition(gene.start, minPos, maxPos, leftPad, usableWidth);
        const x2 = scalePosition(gene.end, minPos, maxPos, leftPad, usableWidth);
        addSvgLine(svg, x1, y, x2, y, "#334155", 2);
        addSvgRect(svg, x1, y - 5, Math.max(4, x2 - x1), 10, gene.strand === "+" ? "#111827" : "#64748b", "none", 0.85);
        addSvgText(svg, x1, y - 8, gene.symbol, "11", "500", "#334155");
      });
    });

    addSvgText(svg, leftPad, mspTop - 8, "MSP regional ancestry composition", "13", "700", "#18212f");
    addSvgText(svg, leftPad, geneTop - 18, "Genes overlapping selected locus", "13", "700", "#18212f");
  }

  function renderContributionPlots() {
    if (!state.qSummary || !window.Plotly) {
      return;
    }
    const globalTrace = {
      type: "bar",
      x: state.qSummary.global.map((entry) => entry.label),
      y: state.qSummary.global.map((entry) => entry.value),
      marker: { color: state.qSummary.global.map((entry) => entry.color) },
      text: state.qSummary.global.map((entry) => entry.value.toFixed(3)),
      textposition: "outside",
      hovertemplate: "%{x}: %{y:.3f}<extra></extra>",
    };
    Plotly.newPlot(ui.globalContributionPlot, [globalTrace], {
      height: 380,
      margin: { l: 40, r: 20, t: 20, b: 40 },
      yaxis: { title: "Mean contribution", range: [0, 1] },
      xaxis: { title: "Ancestry" },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
    }, { responsive: true, displayModeBar: false });

    const locus = selectedLocus();
    const chrom = selectedChromosome();
    const locusSummary = locus
      ? summarizeLocus(locus, mspRowsForChromosome(chrom), genesForChromosome(chrom), state.analysis.ancestries)
      : null;

    const locusTrace = {
      type: "bar",
      x: state.analysis.ancestries.map((entry) => entry.label),
      y: locusSummary ? locusSummary.proportions.map((entry) => entry.value) : state.analysis.ancestries.map(() => 0),
      marker: { color: state.analysis.ancestries.map((entry) => entry.color) },
      text: locusSummary ? locusSummary.proportions.map((entry) => entry.value.toFixed(3)) : state.analysis.ancestries.map(() => "0.000"),
      textposition: "outside",
      hovertemplate: "%{x}: %{y:.3f}<extra></extra>",
    };
    Plotly.newPlot(ui.locusContributionPlot, [locusTrace], {
      height: 380,
      margin: { l: 40, r: 20, t: 36, b: 40 },
      title: {
        text: locusSummary ? `${locusSummary.locus} on chr${locusSummary.chrom}` : "No locus selected",
        font: { size: 14 },
      },
      yaxis: { title: "Mean contribution", range: [0, 1] },
      xaxis: { title: "Ancestry" },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
    }, { responsive: true, displayModeBar: false });
  }

  function renderSamplePlots() {
    if (!state.qSummary || !window.Plotly) {
      return;
    }
    const sortCode = ui.sortAncestrySelect.value || String(state.analysis.ancestries[0].code);
    const samples = state.qSummary.samples
      .slice()
      .sort((a, b) => contributionForCode(b, sortCode) - contributionForCode(a, sortCode));

    const heatmap = {
      type: "heatmap",
      x: state.analysis.ancestries.map((entry) => `${entry.label} (${entry.code})`),
      y: samples.map((entry) => entry.sample),
      z: samples.map((sample) => state.analysis.ancestries.map((entry) => contributionForCode(sample, String(entry.code)))),
      zmin: 0,
      zmax: 1,
      colorscale: [
        [0, "#1a9850"],
        [0.5, "#ffd84d"],
        [1, "#8b0000"],
      ],
      hovertemplate: "Sample=%{y}<br>Ancestry=%{x}<br>Contribution=%{z:.3f}<extra></extra>",
    };
    Plotly.newPlot(ui.sampleHeatmapPlot, [heatmap], {
      height: Math.max(420, Math.min(1200, samples.length * 10)),
      margin: { l: 120, r: 20, t: 20, b: 60 },
      yaxis: { showticklabels: false, ticks: "", title: "Samples" },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
    }, { responsive: true, displayModeBar: false });

    const stackedTraces = state.analysis.ancestries.map((entry) => ({
      type: "bar",
      orientation: "h",
      name: `${entry.label} (${entry.code})`,
      y: samples.map((item) => item.sample),
      x: samples.map((item) => contributionForCode(item, String(entry.code))),
      marker: { color: entry.color },
      hovertemplate: `${entry.label} (code ${entry.code})<br>Sample=%{y}<br>Contribution=%{x:.3f}<extra></extra>`,
    }));
    Plotly.newPlot(ui.sampleStackedPlot, stackedTraces, {
      barmode: "stack",
      height: Math.max(500, Math.min(1400, samples.length * 10)),
      margin: { l: 170, r: 20, t: 20, b: 40 },
      xaxis: { title: "Contribution", range: [0, 1] },
      yaxis: { automargin: false, showticklabels: false, ticks: "", title: "Samples" },
      paper_bgcolor: "#ffffff",
      plot_bgcolor: "#ffffff",
      legend: { orientation: "h", y: 1.08 },
    }, { responsive: true, displayModeBar: false });
  }

  function renderLocusTable() {
    const chrom = selectedChromosome();
    const loci = lociForChromosome(chrom);
    const genes = genesForChromosome(chrom);
    const rows = mspRowsForChromosome(chrom);
    const selected = selectedLocus();

    if (!loci.length) {
      ui.locusTableBody.innerHTML = `<tr><td class="empty" colspan="7">No loci found for chr${escapeHtml(chrom)}.</td></tr>`;
      return;
    }

    const summaryRows = loci.map((locus) => summarizeLocus(locus, rows, genes, state.analysis.ancestries));
    ui.locusTableBody.innerHTML = summaryRows
      .map((row) => {
        const proportions = row.proportions.map((entry) => `${entry.label}:${entry.value.toFixed(3)}`).join(" | ");
        const isSelected = selected && selected.locus === row.locus;
        return `
          <tr${isSelected ? ' style="background:#eef6ff"' : ""}>
            <td>${escapeHtml(row.locus)}</td>
            <td>${escapeHtml(row.chrom)}</td>
            <td>${formatNumber(row.start)}</td>
            <td>${formatNumber(row.end)}</td>
            <td>${escapeHtml(row.genes.join(", ") || "-")}</td>
            <td>${escapeHtml(row.dominant ? `${row.dominant.label} (${row.dominant.code})` : "-")}</td>
            <td>${escapeHtml(proportions)}</td>
          </tr>`;
      })
      .join("");
  }

  function summarizeLocus(locus, mspRows, genes, ancestries) {
    const overlapRows = filterRowsWithinSingleLocus(mspRows, locus);
    const weighted = ancestries.map((entry) => ({ ...entry, sum: 0 }));
    let totalWidth = 0;
    overlapRows.forEach((row) => {
      const overlap = overlapWidth(row, locus);
      if (overlap <= 0) {
        return;
      }
      totalWidth += overlap;
      weighted.forEach((item) => {
        const prop = row.proportions.find((value) => value.code === item.code)?.value || 0;
        item.sum += prop * overlap;
      });
    });
    const proportions = weighted.map((item) => ({
      label: item.label,
      code: item.code,
      color: item.color,
      value: totalWidth ? item.sum / totalWidth : 0,
    }));
    const dominant = proportions.slice().sort((a, b) => b.value - a.value)[0];
    const locusGenes = genes.filter((gene) => overlaps(gene.start, gene.end, locus.start, locus.end)).map((gene) => gene.symbol);
    return {
      ...locus,
      genes: unique(locusGenes),
      dominant,
      proportions,
    };
  }

  function lociForChromosome(chrom) {
    return state.loci.filter((item) => item.chrom === chrom).sort((a, b) => a.start - b.start);
  }

  function genesForChromosome(chrom) {
    return state.genes.filter((item) => item.chrom === chrom).sort((a, b) => a.start - b.start);
  }

  function uniqueGenesWithinLoci(loci, genes) {
    return genes.filter((gene) => loci.some((locus) => overlaps(gene.start, gene.end, locus.start, locus.end)));
  }

  function mspRowsForChromosome(chrom) {
    if (!state.msp) {
      return [];
    }
    return state.msp.rows.filter((item) => item.chrom === chrom).sort((a, b) => a.start - b.start);
  }

  function filterRowsWithinLoci(rows, loci) {
    return rows.filter((row) => loci.some((locus) => overlapWidth(row, locus) > 0));
  }

  function filterRowsWithinSingleLocus(rows, locus) {
    return rows.filter((row) => overlapWidth(row, locus) > 0);
  }

  function clipRowToLoci(row, loci) {
    const segments = [];
    loci.forEach((locus) => {
      const start = Math.max(row.start, locus.start);
      const end = Math.min(row.end, locus.end);
      if (end > start) {
        segments.push({ start, end });
      }
    });
    return segments;
  }

  function layoutGenes(genes) {
    const lanes = [];
    genes.forEach((gene) => {
      let placed = false;
      for (const lane of lanes) {
        const last = lane[lane.length - 1];
        if (last.end < gene.start) {
          lane.push(gene);
          placed = true;
          break;
        }
      }
      if (!placed) {
        lanes.push([gene]);
      }
    });
    return { lanes, laneCount: Math.max(1, lanes.length) };
  }

  function selectedChromosome() {
    return ui.chromosomeSelect.value;
  }

  function selectedLocus() {
    const chrom = selectedChromosome();
    const loci = lociForChromosome(chrom);
    if (!loci.length) {
      return null;
    }
    const selectedName = ui.locusSelect.value;
    return loci.find((entry) => entry.locus === selectedName) || loci[0];
  }

  function updateLocusOptions() {
    const chrom = selectedChromosome();
    const loci = lociForChromosome(chrom);
    fillSelect(
      ui.locusSelect,
      loci.map((locus) => ({
        value: locus.locus,
        label: `${locus.locus} (${formatCompact(locus.start)}-${formatCompact(locus.end)})`,
      })),
      loci[0] ? loci[0].locus : ""
    );
  }

  function contributionForCode(sampleEntry, code) {
    return sampleEntry.contributions.find((item) => String(item.code) === String(code))?.value || 0;
  }

  function tokenize(line) {
    return line.trim().split(/\s+/);
  }

  function requireFile(file, label) {
    if (!file) {
      throw new Error(`Please choose ${label}.`);
    }
    return file;
  }

  function fetchText(path) {
    return fetch(path).then((response) => {
      if (!response.ok) {
        throw new Error(`Failed to fetch ${path}`);
      }
      return response.text();
    });
  }

  function readFileText(file) {
    return new Promise((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(String(reader.result || ""));
      reader.onerror = () => reject(new Error(`Failed to read ${file.name}`));
      reader.readAsText(file);
    });
  }

  function setStatus(message, isError) {
    ui.status.textContent = message;
    ui.status.style.color = isError ? "#b42318" : "#5f6b7a";
  }

  function fillSelect(select, options, selectedValue) {
    select.innerHTML = options
      .map((option) => `<option value="${escapeHtml(option.value)}">${escapeHtml(option.label)}</option>`)
      .join("");
    if (selectedValue != null) {
      select.value = String(selectedValue);
    }
  }

  function clearNode(node) {
    while (node.firstChild) {
      node.removeChild(node.firstChild);
    }
  }

  function drawSvgMessage(svg, message) {
    addSvgText(svg, 40, 60, message, "16", "500", "#5f6b7a");
  }

  function addSvgRect(svg, x, y, width, height, fill, stroke, opacity) {
    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
    rect.setAttribute("x", x);
    rect.setAttribute("y", y);
    rect.setAttribute("width", width);
    rect.setAttribute("height", height);
    rect.setAttribute("fill", fill);
    rect.setAttribute("stroke", stroke || "none");
    if (opacity != null) {
      rect.setAttribute("opacity", opacity);
    }
    svg.appendChild(rect);
  }

  function addSvgLine(svg, x1, y1, x2, y2, stroke, width) {
    const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
    line.setAttribute("x1", x1);
    line.setAttribute("y1", y1);
    line.setAttribute("x2", x2);
    line.setAttribute("y2", y2);
    line.setAttribute("stroke", stroke);
    line.setAttribute("stroke-width", width || 1);
    svg.appendChild(line);
  }

  function addSvgText(svg, x, y, text, size, weight, fill) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", "text");
    node.setAttribute("x", x);
    node.setAttribute("y", y);
    node.setAttribute("font-size", size || "12");
    node.setAttribute("font-weight", weight || "400");
    node.setAttribute("fill", fill || "#111827");
    node.textContent = text;
    svg.appendChild(node);
  }

  function drawAxis(svg, minPos, maxPos, leftPad, usableWidth, y) {
    const tickCount = 6;
    addSvgLine(svg, leftPad, y, leftPad + usableWidth, y, "#64748b", 1);
    for (let i = 0; i <= tickCount; i += 1) {
      const ratio = i / tickCount;
      const value = Math.round(minPos + (maxPos - minPos) * ratio);
      const x = leftPad + usableWidth * ratio;
      addSvgLine(svg, x, y, x, y + 6, "#64748b", 1);
      addSvgText(svg, x - 16, y + 20, formatCompact(value), "11", "400", "#5f6b7a");
    }
  }

  function scalePosition(value, minPos, maxPos, leftPad, usableWidth) {
    if (maxPos === minPos) {
      return leftPad;
    }
    return leftPad + ((value - minPos) / (maxPos - minPos)) * usableWidth;
  }

  function overlapWidth(a, b) {
    return Math.max(0, Math.min(a.end, b.end) - Math.max(a.start, b.start));
  }

  function overlaps(startA, endA, startB, endB) {
    return Math.min(endA, endB) > Math.max(startA, startB);
  }

  function mean(values) {
    return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : 0;
  }

  function chromOrder(chrom) {
    const normalized = String(chrom).replace(/^chr/i, "");
    const numeric = Number(normalized);
    if (!Number.isNaN(numeric)) {
      return numeric;
    }
    return normalized === "X" ? 23 : normalized === "Y" ? 24 : 99;
  }

  function unique(values) {
    return Array.from(new Set(values));
  }

  function formatNumber(value) {
    return new Intl.NumberFormat("en-US").format(value);
  }

  function formatCompact(value) {
    return new Intl.NumberFormat("en-US", { notation: "compact", maximumFractionDigits: 1 }).format(value);
  }

  function escapeHtml(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/\"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }
})();