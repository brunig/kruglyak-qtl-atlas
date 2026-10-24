// Render test for docs/index.html.
//   node tests/render.mjs            (from the repository root; needs Node and Playwright with Chromium)
// Serves docs/ on a local port, screenshots the genome view, a region view and a QTL card in
// light and dark, runs every search type and asserts a result, exercises the table filters and
// TSV download, fails on any console error, and checks a 375 px viewport has no horizontal scroll.
// Google Fonts requests are answered with an empty stylesheet so the test runs offline; the page
// falls back to its system font stack.
import { createRequire } from "node:module";
import { spawn } from "node:child_process";
import { mkdirSync, readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const { chromium } = createRequire(import.meta.url)("playwright");
const ROOT = path.dirname(path.dirname(fileURLToPath(import.meta.url)));
const SHOTS = path.join(ROOT, "tests", "screenshots");
const PORT = 8765 + Math.floor(Math.random() * 500);
const BASE = `http://127.0.0.1:${PORT}/index.html`;
mkdirSync(SHOTS, { recursive: true });

const server = spawn("python3", ["-I", "-m", "http.server", String(PORT), "--bind", "127.0.0.1", "--directory", path.join(ROOT, "docs")], { stdio: "ignore" });
const results = [];
let failed = 0;
function check(name, ok, detail = "") { results.push(`${ok ? "PASS" : "FAIL"}  ${name}${detail ? "  (" + detail + ")" : ""}`); if (!ok) failed++; }

async function waitServer() { for (let i = 0; i < 50; i++) { try { const r = await fetch(BASE); if (r.ok) return; } catch (e) {} await new Promise(r => setTimeout(r, 100)); } throw new Error("server did not start"); }

async function newPage(browser, opts) {
  const ctx = await browser.newContext({ viewport: opts.viewport || { width: 1400, height: 950 }, colorScheme: opts.scheme || "light", acceptDownloads: true, deviceScaleFactor: opts.dpr || 1 });
  await ctx.route(/fonts\.(googleapis|gstatic)\.com/, route => route.fulfill({ status: 200, contentType: "text/css", body: "" }));
  const page = await ctx.newPage();
  page.errors = [];
  page.on("console", m => { if (m.type() === "error") page.errors.push(m.text()); });
  page.on("pageerror", e => page.errors.push(String(e)));
  return { ctx, page };
}
async function load(page, hash = "") {
  await page.goto("about:blank");
  await page.goto(BASE + hash);
  await page.waitForSelector("#loading", { state: "hidden", timeout: 20000 });
  await page.waitForTimeout(150);
}
async function search(page, text) {
  await page.fill("#q", "");
  await page.type("#q", text, { delay: 10 });
  await page.waitForSelector("#sugg:not([hidden]) button", { timeout: 5000 });
  const first = await page.textContent("#sugg button[aria-selected='true']");
  await page.press("#q", "Enter");
  await page.waitForTimeout(250);
  return first;
}

const data = JSON.parse(readFileSync(path.join(ROOT, "docs", "data.json"), "utf8"));
// a within-cross QTL with an assigned candidate gene (PCA1, cadmium) for the card screenshot
const gPCA1 = data.genes.n.indexOf("PCA1");
const qCard = data.qtl.id.findIndex((id, i) => data.qtl.m[i] === "w" && data.qtl.g[i] === gPCA1);
const qJoint = data.qtl.id.findIndex((id, i) => data.qtl.m[i] === "j" && data.qtl.g[i] === gPCA1);

try {
  await waitServer();
  const browser = await chromium.launch();
  for (const scheme of ["light", "dark"]) {
    const { ctx, page } = await newPage(browser, { scheme, dpr: 2 });
    await load(page);
    check(`${scheme}: genome view renders`, await page.isVisible("#gview") && (await page.textContent("#g-title")).includes("Cadmium"));
    await page.screenshot({ path: path.join(SHOTS, `${scheme}-genome.png`), fullPage: true });
    // region view around PCA1 on chrII
    await load(page, "#chrII.760000-830000~t" + data.traits.id[data.traits.raw.indexOf("Cadmium_Chloride;75uM;2")]);
    const rows = await page.$$eval("#ltable tr.row", r => r.length);
    check(`${scheme}: region view lists QTL`, rows > 0, rows + " rows");
    check(`${scheme}: region shows candidate gene PCA1`, (await page.textContent("#candlist")).includes("PCA1"));
    await page.screenshot({ path: path.join(SHOTS, `${scheme}-region.png`), fullPage: true });
    // QTL card with its source row
    await load(page, "#chrII.760000-830000~q" + data.qtl.id[qCard]);
    await page.waitForSelector("#qcard:not([hidden])");
    await page.click("#qc-src");
    await page.waitForSelector("#qc-srcbox table", { timeout: 5000 });
    const card = await page.textContent("#qcard");
    check(`${scheme}: QTL card shows trait, LOD, effect, candidate, source row`, card.includes("Cadmium") && card.includes("LOD") && card.includes("Increased by") && card.includes("PCA1") && card.includes("within-cross model"));
    await page.addStyleTag({ content: ".bar { position: static !important; }" });
    await (await page.$("#qcard")).screenshot({ path: path.join(SHOTS, `${scheme}-card.png`) });
    // a joint-model card
    await load(page, "#chrII.760000-830000~q" + data.qtl.id[qJoint]);
    const jcard = await page.textContent("#qcard");
    check(`${scheme}: joint QTL card shows per-cross effects and allele`, jcard.includes("Effect in each cross") && /alternate allele|reference \(S288C\) allele/.test(jcard));
    await page.addStyleTag({ content: ".bar { position: static !important; }" });
    if (scheme === "light") await (await page.$("#qcard")).screenshot({ path: path.join(SHOTS, `light-card-joint.png`) });
    // table view
    await load(page, "#table");
    await page.screenshot({ path: path.join(SHOTS, `${scheme}-table.png`) });
    check(`${scheme}: no console errors`, page.errors.length === 0, page.errors.join(" | "));
    await ctx.close();
  }

  // ---- searches and table, light desktop ----
  {
    const { ctx, page } = await newPage(browser, {});
    await load(page);
    let f = await search(page, "IRA2");
    check("search gene name IRA2", f.includes("IRA2") && await page.isVisible("#lview") && (await page.textContent("#genecard")).includes("IRA2"));
    f = await search(page, "YOL081W");
    check("search systematic name YOL081W", f.includes("IRA2") && (await page.textContent("#genecard")).includes("YOL081W"));
    f = await search(page, "chrIV:509600");
    check("search coordinate chrIV:509600", (await page.textContent("#l-title")) === "chrIV:459,600–559,600");
    f = await search(page, "IV:400k-600k");
    check("search range IV:400k-600k", (await page.textContent("#l-title")) === "chrIV:400,000–600,000");
    f = await search(page, "caffeine");
    check("search trait caffeine", f.includes("Caffeine") && await page.isVisible("#gview") && (await page.textContent("#g-title")).includes("Caffeine"));
    f = await search(page, "375");
    check("search cross code 375", f.includes("M22 × BY") && (await page.textContent("#crossf")).includes("M22 × BY"));
    f = await search(page, "M22 x BY");
    check("search cross by parents 'M22 x BY'", f.includes("M22 × BY"));
    f = await search(page, "YJM454");
    check("search parent strain YJM454", f.includes("YJM454") && (await page.textContent("#crossf")).includes("2 crosses"));
    f = await search(page, "round robin");
    check("search dataset 'round robin'", f.includes("Bloom 2019") && (await page.getAttribute("#dschips .chip", "aria-pressed")) === "true");
    f = await search(page, "Bloom");
    check("search dataset 'Bloom'", f.includes("Bloom 2019"));
    // dataset chip toggling (single dataset cannot be switched off)
    await page.click("#dschips .chip");
    check("dataset chip stays on when it is the only one", (await page.getAttribute("#dschips .chip", "aria-pressed")) === "true");
    // model control
    await load(page, "#genome");
    await page.click("#m-j"); await page.waitForTimeout(150);
    check("model control: joint only", (await page.textContent("#g-meta")).includes(" 0 within-cross"));
    await page.click("#m-both");
    // table: filters, sort, download
    await load(page, "#table");
    const total = await page.textContent("#tv-title");
    check("table shows all 12,303 QTL", total.includes("12,303 of 12,303"), total);
    const caf = data.traits.raw.indexOf("Caffeine;15mM;2");
    await page.selectOption("#tf-t", String(caf)); await page.waitForTimeout(150);
    await page.selectOption("#tf-x", "-1"); await page.waitForTimeout(150);
    const t2 = await page.textContent("#tv-title");
    const nCafJ = data.qtl.t.filter((t, i) => t === caf && data.qtl.m[i] === "j").length;
    check("table filter trait + joint", t2.includes(`${nCafJ.toLocaleString("en-US")} of`), t2);
    await page.selectOption("#tf-x", ""); await page.fill("#tf-lod", "20"); await page.waitForTimeout(400);
    const nCaf20 = data.qtl.t.filter((t, i) => t === caf && data.qtl.lod[i] != null && data.qtl.lod[i] >= 20).length;
    check("table filter LOD >= 20", (await page.textContent("#tv-title")).includes(`${nCaf20.toLocaleString("en-US")} of`));
    await page.click("#ttable th[data-sort='lod']");
    const firstLod = await page.$eval("#ttable tbody tr td:nth-child(10)", td => parseFloat(td.textContent));
    const maxLod = Math.max(...data.qtl.lod.filter((l, i) => l != null && data.qtl.t[i] === caf));
    check("table sorts by LOD", Math.abs(firstLod - maxLod) < 0.06, firstLod + " vs " + maxLod);
    const [dl] = await Promise.all([page.waitForEvent("download"), page.click("#tsv")]);
    const tsv = readFileSync(await dl.path(), "utf8").trim().split("\n");
    check("TSV download has header + filtered rows", tsv.length === nCaf20 + 1 && tsv[0].startsWith("id\tdataset\ttrait"), (tsv.length - 1) + " rows");
    // hash state round trip
    await load(page, "#chrXV.500000-620000~mw~x" + data.crosses.id[0]);
    check("hash restores region, model and cross", (await page.textContent("#l-title")).startsWith("chrXV:500,000") && (await page.getAttribute("#m-w", "aria-pressed")) === "true" && (await page.textContent("#crossf")).includes("M22 × BY"));
    check("no console errors (searches, table)", page.errors.length === 0, page.errors.join(" | "));
    await ctx.close();
  }

  // ---- phone width ----
  {
    const { ctx, page } = await newPage(browser, { viewport: { width: 375, height: 812 }, dpr: 2 });
    for (const [name, hash] of [["genome", "#genome"], ["region", "#chrII.760000-830000~q" + data.qtl.id[qCard]], ["table", "#table"]]) {
      await load(page, hash);
      const sw = await page.evaluate(() => Math.max(document.documentElement.scrollWidth, document.body.scrollWidth));
      check(`375 px ${name}: no horizontal scroll`, sw <= 375, "scrollWidth " + sw);
      await page.screenshot({ path: path.join(SHOTS, `phone-${name}.png`), fullPage: name !== "table" });
    }
    check("375 px: no console errors", page.errors.length === 0, page.errors.join(" | "));
    await ctx.close();
  }
  await browser.close();
} catch (e) {
  check("test harness", false, String(e && e.stack || e));
} finally {
  server.kill();
}
console.log(results.join("\n"));
console.log(`\n${results.length - failed} passed, ${failed} failed`);
process.exit(failed ? 1 : 0);
