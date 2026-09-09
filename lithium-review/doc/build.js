const fs = require("fs");
const path = require("path");
const D = require("docx");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle,
  ImageRun, PageBreak, PageNumber, Footer, Header, TableOfContents,
  PageOrientation, convertInchesToTwip, LevelFormat, TabStopType,
} = D;

const HERE = __dirname;
const FIGDIR = path.join(HERE, "..", "analysis", "figures");
const CSV = path.join(HERE, "..", "analysis", "li_recovery_dataset.csv");

const TABLES = require("./tables.js");
const REFS = require("./refs.js");
const BODY = [].concat(require("./body1.js"), require("./body2.js"));

// ---------------------------------------------------------------- constants
const PAGE_W = 12240, PAGE_H = 15840;      // US Letter, DXA
const MARGIN = convertInchesToTwip(1);
const CONTENT_W = PAGE_W - 2 * MARGIN;      // 9360
const SERIF = "Times New Roman";
const GREY = "F0F0F0", RULE = "9A9A9A";

let figureCount = 0, tableCount = 0;

// ---------------------------------------------------------------- helpers
function body(text, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.JUSTIFIED,
    spacing: { after: opts.after ?? 140, line: opts.line ?? 300 },
    indent: opts.indent,
    children: [new TextRun({ text, font: SERIF, size: opts.size || 21,
      italics: !!opts.italics, bold: !!opts.bold, color: opts.color })],
  });
}

function h(text, level) {
  const sizes = { 1: 28, 2: 24, 3: 22 };
  return new Paragraph({
    heading: level === 1 ? HeadingLevel.HEADING_1
           : level === 2 ? HeadingLevel.HEADING_2 : HeadingLevel.HEADING_3,
    spacing: { before: level === 1 ? 360 : 260, after: level === 1 ? 180 : 120 },
    children: [new TextRun({ text, font: SERIF, size: sizes[level], bold: true, color: "1A1A1A" })],
  });
}

function bullet(text) {
  return new Paragraph({
    numbering: { reference: "bullets", level: 0 },
    alignment: AlignmentType.JUSTIFIED,
    spacing: { after: 90, line: 290 },
    children: [new TextRun({ text, font: SERIF, size: 21 })],
  });
}

function caption(text) {
  return new Paragraph({
    alignment: AlignmentType.LEFT,
    spacing: { before: 90, after: 240, line: 260 },
    children: [new TextRun({ text, font: SERIF, size: 18, italics: true, color: "3A3A3A" })],
  });
}

function cellPara(text, { bold = false, align = AlignmentType.LEFT, size = 17 } = {}) {
  return new Paragraph({
    alignment: align,
    spacing: { before: 40, after: 40, line: 240 },
    children: [new TextRun({ text: String(text), font: SERIF, size, bold })],
  });
}

// scale a width array so it exactly fills the content width
function scale(widths, total = CONTENT_W) {
  const sum = widths.reduce((a, b) => a + b, 0);
  const out = widths.map(w => Math.floor(w * total / sum));
  out[out.length - 1] += total - out.reduce((a, b) => a + b, 0);
  return out;
}

function makeTable(spec) {
  tableCount += 1;
  const widths = scale(spec.widths);
  const border = { style: BorderStyle.SINGLE, size: 2, color: RULE };
  const noneB = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };

  const headRow = new TableRow({
    tableHeader: true,
    children: spec.head.map((txt, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: { type: ShadingType.CLEAR, fill: GREY, color: "auto" },
      margins: { top: 60, bottom: 60, left: 90, right: 90 },
      borders: { top: border, bottom: border, left: noneB, right: noneB },
      children: [cellPara(txt, { bold: true, size: 17 })],
    })),
  });

  const rows = spec.rows.map((r, ri) => new TableRow({
    cantSplit: true,
    children: r.map((txt, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: ri % 2 === 1
        ? { type: ShadingType.CLEAR, fill: "F9F9F9", color: "auto" }
        : undefined,
      margins: { top: 50, bottom: 50, left: 90, right: 90 },
      borders: {
        top: noneB,
        bottom: ri === spec.rows.length - 1 ? border : { style: BorderStyle.SINGLE, size: 1, color: "DDDDDD" },
        left: noneB, right: noneB,
      },
      children: [cellPara(txt, { align: i === 0 ? AlignmentType.LEFT : AlignmentType.LEFT })],
    })),
  }));

  const out = [
    new Paragraph({
      spacing: { before: 260, after: 100 },
      children: [new TextRun({ text: spec.title, font: SERIF, size: 19, bold: true })],
    }),
    new Table({
      columnWidths: widths,
      width: { size: CONTENT_W, type: WidthType.DXA },
      rows: [headRow, ...rows],
    }),
  ];
  if (spec.note) {
    out.push(new Paragraph({
      spacing: { before: 90, after: 260, line: 250 },
      children: [new TextRun({ text: spec.note, font: SERIF, size: 16, italics: true, color: "555555" })],
    }));
  } else {
    out.push(new Paragraph({ spacing: { after: 200 }, children: [] }));
  }
  return out;
}

function makeFigure(spec) {
  figureCount += 1;
  const file = path.join(FIGDIR, spec.file);
  if (!fs.existsSync(file)) {
    console.warn("missing figure:", spec.file);
    return [body("[figure unavailable: " + spec.file + "]", { italics: true })];
  }
  // native px at 200 dpi -> scale to fit content width (6.5 in)
  const maxW = 600;          // points-ish; docx uses px at 96dpi for width/height
  const dims = { width: maxW, height: Math.round(maxW * 0.52) };
  // figure-specific aspect ratios
  const ar = { "fig1_capacity.png": 0.545, "fig2_energy.png": 0.567,
               "fig3_recovery.png": 0.531, "fig4_energy_trend.png": 0.567,
               "fig5_selectivity.png": 0.567, "fig6_completeness.png": 0.533,
               "fig7_cycling.png": 0.567 };
  if (ar[spec.file]) dims.height = Math.round(maxW * ar[spec.file]);

  return [
    new Paragraph({
      alignment: AlignmentType.CENTER,
      spacing: { before: 220, after: 40 },
      children: [new ImageRun({
        type: "png",
        data: fs.readFileSync(file),
        transformation: dims,
      })],
    }),
    caption(spec.caption),
  ];
}

// ---------------------------------------------------------------- title page
const titlePage = [
  new Paragraph({ spacing: { before: 1400, after: 0 }, alignment: AlignmentType.CENTER,
    children: [new TextRun({ text: "REVIEW AND META-ANALYSIS", font: SERIF, size: 20,
      bold: true, color: "666666", characterSpacing: 60 })] }),
  new Paragraph({ spacing: { before: 320, after: 200 }, alignment: AlignmentType.CENTER,
    children: [new TextRun({
      text: "Lithium Recovery from Aqueous Resources: A Systematic Review and Quantitative Meta-Analysis of Adsorption, Electrochemical, Membrane, Solvent Extraction and Precipitation Technologies",
      font: SERIF, size: 34, bold: true })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 200, after: 60 },
    children: [new TextRun({ text: "A pooled synthesis of 131 quantitative performance records drawn from six full text reviews spanning eight technology families", font: SERIF, size: 22, italics: true, color: "444444" })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 700, after: 40 },
    children: [new TextRun({ text: "Prepared for M. Al-Shebli", font: SERIF, size: 22 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 40 },
    children: [new TextRun({ text: "9 September 2026", font: SERIF, size: 22 })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 900, after: 60 },
    children: [new TextRun({ text: "Evidence base", font: SERIF, size: 19, bold: true, color: "555555" })] }),
  new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 30, line: 280 },
    children: [new TextRun({ text: "6 full text review articles  |  131 harmonised performance records  |  8 technology families  |  283 references", font: SERIF, size: 19, color: "555555" })] }),
  new Paragraph({ children: [new PageBreak()] }),
];

// ---------------------------------------------------------------- contents
const toc = [
  h("Contents", 1),
  new TableOfContents("Contents", { hyperlink: true, headingStyleRange: "1-2" }),
  new Paragraph({ spacing: { before: 200 },
    children: [new TextRun({ text: "If the entries below appear blank, open the document in Word and press Ctrl+A then F9 to build the field.", font: SERIF, size: 17, italics: true, color: "666666" })] }),
  new Paragraph({ children: [new PageBreak()] }),
];

// ---------------------------------------------------------------- body
const content = [];
for (const blk of BODY) {
  if (blk.h1) content.push(h(blk.h1, 1));
  else if (blk.h2) content.push(h(blk.h2, 2));
  else if (blk.h3) content.push(h(blk.h3, 3));
  else if (blk.p) content.push(body(blk.p));
  else if (blk.bullets) blk.bullets.forEach(b => content.push(bullet(b)));
  else if (blk.table) {
    const spec = TABLES[blk.table];
    if (!spec) { console.warn("missing table:", blk.table); continue; }
    content.push(...makeTable(spec));
  }
  else if (blk.figure) content.push(...makeFigure(blk.figure));
  else if (blk.pagebreak) content.push(new Paragraph({ children: [new PageBreak()] }));
}

// ---------------------------------------------------------------- references
const refBlocks = [
  new Paragraph({ children: [new PageBreak()] }),
  h("References", 1),
  body("References are listed alphabetically within thematic groupings that follow the order in which the material appears in the text. Entries transcribed from the reference lists of the six source reviews retain the bibliographic detail given there. Entries identified by literature search during preparation of this manuscript and whose full text could not be retrieved are given with the identifier that was verified, and volume or page numbers are omitted rather than inferred.", { italics: true, size: 19 }),
];
REFS.forEach((r, i) => {
  refBlocks.push(new Paragraph({
    alignment: AlignmentType.JUSTIFIED,
    spacing: { after: 70, line: 250 },
    indent: { left: 480, hanging: 480 },
    children: [
      new TextRun({ text: `[${i + 1}]  `, font: SERIF, size: 18, bold: true }),
      new TextRun({ text: r, font: SERIF, size: 18 }),
    ],
  }));
});

// ---------------------------------------------------------------- appendix
function parseCSV(txt) {
  const lines = txt.trim().split(/\r?\n/);
  const parse = (line) => {
    const out = []; let cur = "", q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (c === '"') { if (q && line[i + 1] === '"') { cur += '"'; i++; } else q = !q; }
      else if (c === "," && !q) { out.push(cur); cur = ""; }
      else cur += c;
    }
    out.push(cur); return out;
  };
  const head = parse(lines[0]);
  return lines.slice(1).map(l => {
    const c = parse(l); const o = {};
    head.forEach((hh, i) => o[hh] = c[i] ?? "");
    return o;
  });
}

const rows = parseCSV(fs.readFileSync(CSV, "utf8"));
const appendixSpec = {
  title: "Appendix Table A1. Complete harmonised evidence table (131 records).",
  note: "Cap = uptake or extraction capacity (mg g-1); E = specific energy (Wh mol-1 Li); Rec = recovery or extraction efficiency (%); Pur = product purity (%); SF = Li/Mg separation factor; Cyc = cycles tested; Ret = capacity retention (%). Src identifies the source review: R1 Liu et al. (2019); R2 Wang et al. (2022); R3 Zavahir et al. (2021); R4 Wu et al. (2022); R5 Zhao et al. (2019); R6 Khalil et al. (2022). A dash indicates the metric was not reported in the source.",
  head: ["ID", "Study", "Yr", "Family", "Material or system", "Cap", "E", "Rec", "Pur", "SF", "Cyc", "Ret", "Src"],
  widths: [420, 1250, 400, 1450, 2250, 420, 430, 420, 400, 470, 380, 380, 330],
  rows: rows.map(r => [
    r.id,
    r.study.replace(" et al.", " et al.").replace("compiled entry (R3 Table 1)", "compiled, R3 T1").replace("compiled entry (R3 Table 2)", "compiled, R3 T2"),
    r.year || "-",
    r.family.replace("Electrodialysis / membrane electro-process", "Electrodialysis")
            .replace("Electrosorption / ion pumping", "Electrosorption")
            .replace("Adsorption / ion sieve", "Adsorption")
            .replace("Nanofiltration / selective membrane", "Nanofiltration")
            .replace("Hybrid membrane-adsorbent", "Hybrid membrane")
            .replace("Capacitive deionisation", "CDI"),
    (r.material_or_system || "-").slice(0, 46),
    r.capacity_mg_g || "-",
    r.energy_Wh_mol || "-",
    r.recovery_pct || "-",
    r.purity_pct || "-",
    r.sel_Li_Mg || "-",
    r.cycles || "-",
    r.retention_pct || "-",
    r.compiled_in,
  ]),
};

const appendix = [
  new Paragraph({ children: [new PageBreak()] }),
  h("Appendix A. Complete harmonised evidence table", 1),
  body("The table below reproduces every record underlying the pooled statistics reported in Section 5. Each row records the primary study as named in the source review, the technology family to which it was assigned, the material or system tested, and the six target metrics where reported. The final column identifies the source review from which the record was transcribed, so that any value can be traced back to the document it came from and, through that document's reference list, to the primary publication."),
];

// ---------------------------------------------------------------- document
const doc = new Document({
  creator: "Lithium recovery evidence synthesis",
  title: "Lithium Recovery from Aqueous Resources: A Systematic Review and Quantitative Meta-Analysis",
  description: "Pooled meta-analysis of 131 performance records across eight lithium recovery technology families",
  numbering: {
    config: [{
      reference: "bullets",
      levels: [{
        level: 0, format: LevelFormat.BULLET, text: "•",
        alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 520, hanging: 260 } } },
      }],
    }],
  },
  styles: {
    default: {
      document: { run: { font: SERIF, size: 21 } },
      heading1: { run: { font: SERIF, size: 28, bold: true, color: "1A1A1A" } },
      heading2: { run: { font: SERIF, size: 24, bold: true, color: "1A1A1A" } },
      heading3: { run: { font: SERIF, size: 22, bold: true, color: "1A1A1A" } },
    },
  },
  sections: [
    {
      properties: {
        page: { size: { width: PAGE_W, height: PAGE_H },
                margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } },
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [PageNumber.CURRENT], font: SERIF, size: 18, color: "666666" })],
          })],
        }),
      },
      children: [...titlePage, ...toc, ...content, ...refBlocks],
    },
    {
      properties: {
        page: {
          size: { width: PAGE_W, height: PAGE_H, orientation: PageOrientation.LANDSCAPE },
          margin: { top: 720, bottom: 720, left: 720, right: 720 },
        },
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [new TextRun({ children: [PageNumber.CURRENT], font: SERIF, size: 18, color: "666666" })],
          })],
        }),
      },
      children: (() => {
        // landscape section: recompute widths against landscape content width
        const LW = PAGE_H - 1440;   // 15840 - 2*720
        const saved = CONTENT_W;
        const widths = scale(appendixSpec.widths, LW);
        const border = { style: BorderStyle.SINGLE, size: 2, color: RULE };
        const noneB = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
        const headRow = new TableRow({
          tableHeader: true,
          children: appendixSpec.head.map((t, i) => new TableCell({
            width: { size: widths[i], type: WidthType.DXA },
            shading: { type: ShadingType.CLEAR, fill: GREY, color: "auto" },
            margins: { top: 50, bottom: 50, left: 70, right: 70 },
            borders: { top: border, bottom: border, left: noneB, right: noneB },
            children: [cellPara(t, { bold: true, size: 15 })],
          })),
        });
        const trows = appendixSpec.rows.map((r, ri) => new TableRow({
          cantSplit: true,
          children: r.map((t, i) => new TableCell({
            width: { size: widths[i], type: WidthType.DXA },
            shading: ri % 2 === 1 ? { type: ShadingType.CLEAR, fill: "F9F9F9", color: "auto" } : undefined,
            margins: { top: 35, bottom: 35, left: 70, right: 70 },
            borders: { top: noneB, bottom: { style: BorderStyle.SINGLE, size: 1, color: "E2E2E2" },
                       left: noneB, right: noneB },
            children: [cellPara(t, { size: 14 })],
          })),
        }));
        return [
          ...appendix,
          new Paragraph({
            spacing: { before: 160, after: 100 },
            children: [new TextRun({ text: appendixSpec.title, font: SERIF, size: 19, bold: true })],
          }),
          new Table({
            columnWidths: widths,
            width: { size: LW, type: WidthType.DXA },
            rows: [headRow, ...trows],
          }),
          new Paragraph({
            spacing: { before: 100 },
            children: [new TextRun({ text: appendixSpec.note, font: SERIF, size: 15, italics: true, color: "555555" })],
          }),
        ];
      })(),
    },
  ],
});

const OUT = path.join(HERE, "Lithium_Recovery_Review_and_Meta_Analysis.docx");
Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(OUT, buf);
  console.log("written:", OUT, (buf.length / 1024).toFixed(0) + " KB");
  console.log("tables:", tableCount, "figures:", figureCount, "refs:", REFS.length,
              "appendix rows:", rows.length);
});
