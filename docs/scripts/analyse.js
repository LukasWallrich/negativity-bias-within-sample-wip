const XLSX = require('xlsx');
const fs = require('fs');

const root = '/home/lukas/.openclaw/workspace/negativity-bias-wip';
const sheet = XLSX.readFile(`${root}/data/F2.xlsx`).Sheets['Complete raw data set'];
const rows = XLSX.utils.sheet_to_json(sheet, {range: 1, defval: null});

const num = x => typeof x === 'number' ? x : Number(x);
const fisher = r => Math.atanh(Math.max(-0.999999, Math.min(0.999999, r)));
const keyOf = x => [x.Paper_No, x.Study_No, x.Sample_No].join('|');
const clean = x => String(x ?? '').replace(/\s+/g, ' ').trim();

// The workbook's v74 flag marks rows retained for the authors' within-subjects
// screening. We additionally require a shared sample and shared DV code.
const eligible = rows.filter(x => String(x.v74_Within_Subjects_Inc_Discard) === '1'
  && ['0', '1'].includes(String(x.v45_Recode_Pos_Ambiv_Neg_Unc))
  && String(x.v45d_ES_Congruence) === '0'
  && Number.isFinite(num(x.v72a_ES_Per_P_T))
  && Number.isFinite(num(x.v55_Sample_N))
  && x.v49_DV_Details != null);

const bySample = new Map();
for (const x of eligible) {
  const k = keyOf(x);
  if (!bySample.has(k)) bySample.set(k, []);
  bySample.get(k).push(x);
}

function bestPair(sampleRows) {
  const byDv = new Map();
  for (const x of sampleRows) {
    const k = String(x.v49_DV_Details);
    if (!byDv.has(k)) byDv.set(k, {pos: [], neg: []});
    byDv.get(k)[String(x.v45_Recode_Pos_Ambiv_Neg_Unc) === '0' ? 'pos' : 'neg'].push(x);
  }
  const pairs = [];
  for (const [dv, z] of byDv) {
    if (!z.pos.length || !z.neg.length) continue;
    const pos = [...z.pos].sort((a,b) => num(a.v57_Purity_Ax||99)-num(b.v57_Purity_Ax||99))[0];
    const neg = [...z.neg].sort((a,b) => num(a.v57_Purity_Ax||99)-num(b.v57_Purity_Ax||99))[0];
    pairs.push({dv, pos, neg, purity: num(pos.v57_Purity_Ax||99)+num(neg.v57_Purity_Ax||99)});
  }
  return pairs.sort((a,b) => a.purity-b.purity)[0] || null;
}

const pairs = [];
for (const [sample, sampleRows] of bySample) {
  const p = bestPair(sampleRows);
  if (!p) continue;
  const n = Math.min(num(p.pos.v55_Sample_N), num(p.neg.v55_Sample_N));
  const posR = num(p.pos.v72a_ES_Per_P_T);
  const negR = num(p.neg.v72a_ES_Per_P_T);
  // P&T scoring: positive contact effects are negative when beneficial;
  // negative contact effects are positive when detrimental. Re-express both
  // as congruent magnitudes before differencing.
  const zPos = fisher(-posR);
  const zNeg = fisher(negR);
  pairs.push({
    sample, paper: p.pos.Paper_No, study: p.pos.Study_No, sampleNo: p.pos.Sample_No,
    n, dv: p.dv, target: clean(p.pos.v50_Target_OG), ref: clean(p.pos.BibRef),
    pos_r: posR, neg_r: negR, pos_z: zPos, neg_z: zNeg,
    contrast: zNeg-zPos, purity: p.purity,
    pos_desc: clean(p.pos.v47_IV_Details), neg_desc: clean(p.neg.v47_IV_Details)
  });
}

function pool(data, rho) {
  const vi = data.map(x => {
    const a = 1/(x.n-3), b = 1/(x.n-3);
    return a+b-2*rho*Math.sqrt(a*b);
  });
  const wi = vi.map(v => 1/v);
  const fixed = data.reduce((s,x,i) => s+wi[i]*x.contrast,0)/wi.reduce((s,x)=>s+x,0);
  const q = data.reduce((s,x,i) => s+wi[i]*(x.contrast-fixed)**2,0);
  const df = data.length-1;
  const c = wi.reduce((s,x)=>s+x,0)-wi.reduce((s,x)=>s+x*x,0)/wi.reduce((s,x)=>s+x,0);
  const tau2 = Math.max(0,(q-df)/Math.max(c,1e-12));
  const wr = vi.map(v => 1/(v+tau2));
  const est = data.reduce((s,x,i) => s+wr[i]*x.contrast,0)/wr.reduce((s,x)=>s+x,0);
  const se = Math.sqrt(1/wr.reduce((s,x)=>s+x,0));
  return {rho, k:data.length, estimate:est, se, lo:est-1.96*se, hi:est+1.96*se,
    tau2, I2:q>df ? Math.max(0,(q-df)/q) : 0, z:est/se, p:2*(1-normalCdf(Math.abs(est/se))) };
}
function normalCdf(x) { const t=1/(1+0.2316419*x), d=0.3989423*Math.exp(-x*x/2); return 1-d*t*(0.3193815+t*(-0.3565638+t*(1.781478+t*(-1.821256+t*1.330274*t)))); }

const results = [0,.25,.5,.75,.9].map(rho => pool(pairs,rho));
const summary = {
  source_rows: rows.length, flagged_rows: rows.filter(x=>String(x.v74_Within_Subjects_Inc_Discard)==='1').length,
  eligible_rows: eligible.length, paired_samples: pairs.length,
  paired_sample_n: pairs.reduce((s,x)=>s+x.n,0), results,
  caveat: 'Sampling covariance between positive- and negative-contact correlations is not reported; results are sensitivity analyses over assumed covariance rho.'
};
fs.writeFileSync(`${root}/report/pairs.json`, JSON.stringify(pairs,null,2));
fs.writeFileSync(`${root}/report/summary.json`, JSON.stringify(summary,null,2));
console.log(JSON.stringify(summary,null,2));
