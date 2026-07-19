const XLSX = require('xlsx');
const fs = require('fs');

const root = '/home/lukas/.openclaw/workspace/negativity-bias-wip';
const rows = XLSX.utils.sheet_to_json(
  XLSX.readFile(`${root}/data/F2.xlsx`).Sheets['Complete raw data set'],
  {range: 1, defval: null}
);
const num = x => typeof x === 'number' ? x : Number(x);
const clean = x => String(x ?? '').replace(/\s+/g, ' ').trim().toLowerCase();
const fisher = r => Math.atanh(Math.max(-.999999, Math.min(.999999, r)));
const sampleKey = x => [x.Paper_No, x.Study_No, x.Sample_No].join('|');
const outcomeKey = x => `${sampleKey(x)}||${clean(x.v47_IV_Details)}`;

// Keep the same conceptual screen as the published exploratory report, but
// retain every eligible effect size in a source-coded outcome group that has
// both positive and negative contact. This avoids choosing one "purest" pair.
const eligible = rows.filter(x =>
  String(x.v74_Within_Subjects_Inc_Discard) === '1' &&
  ['0', '1'].includes(String(x.v45_Recode_Pos_Ambiv_Neg_Unc)) &&
  String(x.v45d_ES_Congruence) === '0' &&
  Number.isFinite(num(x.v72a_ES_Per_P_T)) &&
  Number.isFinite(num(x.v55_Sample_N)) &&
  x.v49_DV_Details != null
);
const grouped = new Map();
for (const x of eligible) {
  const k = outcomeKey(x);
  if (!grouped.has(k)) grouped.set(k, []);
  grouped.get(k).push(x);
}
const matchedGroups = [...grouped.values()].filter(g =>
  g.some(x => String(x.v45_Recode_Pos_Ambiv_Neg_Unc) === '0') &&
  g.some(x => String(x.v45_Recode_Pos_Ambiv_Neg_Unc) === '1')
);
const all = matchedGroups.flat().map(x => ({
  y: fisher(num(x.v72a_ES_Per_P_T)),
  vi: 1 / Math.max(1, num(x.v55_Sample_N) - 3),
  neg: String(x.v45_Recode_Pos_Ambiv_Neg_Unc) === '1' ? 1 : 0,
  highPrej: (x.v8_Prej_Raw === 1 || x.v8_Prej_Raw === 3) ? 1 : 0,
  classifiedPrej: [1, 2, 3].includes(num(x.v8_Prej_Raw)),
  sample: sampleKey(x), paper: x.Paper_No, study: x.Study_No,
  raw: x
}));

function normalCdf(x) {
  const a = Math.abs(x), t = 1/(1+0.2316419*a), d=.3989423*Math.exp(-a*a/2);
  const p=1-d*t*(.3193815+t*(-.3565638+t*(1.781478+t*(-1.821256+t*1.330274*t))));
  return x >= 0 ? p : 1-p;
}
function inv(A) {
  const n=A.length, M=A.map((r,i)=>r.slice().concat(Array.from({length:n},(_,j)=>i===j?1:0)));
  for(let c=0;c<n;c++) { let p=c; for(let i=c+1;i<n;i++) if(Math.abs(M[i][c])>Math.abs(M[p][c])) p=i;
    if(Math.abs(M[p][c])<1e-12) throw new Error('Singular matrix'); [M[c],M[p]]=[M[p],M[c]];
    const d=M[c][c]; for(let j=0;j<2*n;j++) M[c][j]/=d;
    for(let i=0;i<n;i++) if(i!==c){const f=M[i][c]; for(let j=0;j<2*n;j++) M[i][j]-=f*M[c][j];}
  } return M.map(r=>r.slice(n));
}
function dot(a,b){return a.reduce((s,x,i)=>s+x*b[i],0)}
function wls(data, columns, label) {
  const p=columns[0].length, XtWX=Array.from({length:p},()=>Array(p).fill(0)), XtWy=Array(p).fill(0), groups=new Map();
  data.forEach((d,i)=>{const x=columns[i], w=1/d.vi; for(let a=0;a<p;a++){XtWy[a]+=w*x[a]*d.y;for(let b=0;b<p;b++)XtWX[a][b]+=w*x[a]*x[b];} if(!groups.has(d.sample))groups.set(d.sample,[]);groups.get(d.sample).push({x,y:d.y,w});});
  const V=inv(XtWX), beta=V.map(r=>dot(r,XtWy));
  const meat=Array.from({length:p},()=>Array(p).fill(0));
  for(const g of groups.values()){const s=Array(p).fill(0);for(const q of g){const e=q.y-dot(q.x,beta);for(let a=0;a<p;a++)s[a]+=q.w*q.x[a]*e;}for(let a=0;a<p;a++)for(let b=0;b<p;b++)meat[a][b]+=s[a]*s[b];}
  const cov=Array.from({length:p},()=>Array(p).fill(0));
  for(let a=0;a<p;a++) for(let b=0;b<p;b++)
    for(let k=0;k<p;k++) for(let l=0;l<p;l++) cov[a][b]+=V[a][k]*meat[k][l]*V[l][b];
  const terms=beta.map((b,i)=>({estimate:b,se:Math.sqrt(Math.max(0,cov[i][i])),z:b/Math.sqrt(Math.max(1e-30,cov[i][i]))}));
  terms.forEach(t=>{t.lo=t.estimate-1.96*t.se;t.hi=t.estimate+1.96*t.se;t.p=2*(1-normalCdf(Math.abs(t.z)))});
  return {label,n:data.length,samples:groups.size,terms,cov};
}

const colsBasic=all.map(d=>[1,d.neg]);
const basic=wls(all,colsBasic,'all rows: z(r) ~ negative contact');
const basicAsymmetry=linearCombo(basic,[2,1],'negative minus reversed positive; all rows');
const classified=all.filter(d=>d.classifiedPrej);
const colsPrej=classified.map(d=>[1,d.neg,d.highPrej,d.neg*d.highPrej]);
const prejudice=wls(classified,colsPrej,'classified prejudice: z(r) ~ negative * high prejudice');

function linearCombo(model, weights, label) {
  const estimate=dot(weights,model.terms.map(t=>t.estimate));
  const variance=weights.reduce((s,_,i)=>s+weights.reduce((t,__,j)=>t+weights[i]*model.cov[i][j]*weights[j],0),0);
  const se=Math.sqrt(Math.max(0,variance)), z=estimate/Math.max(se,1e-30);
  return {label,estimate,se,lo:estimate-1.96*se,hi:estimate+1.96*se,z,p:2*(1-normalCdf(Math.abs(z)))};
}
// The moderator coefficient itself is not the asymmetry contrast because the
// positive-contact estimate is expected to be negative. These are the paired
// contrasts implied by the all-row model: z(negative) - z(-positive).
const asymmetry={non_prejudiced:linearCombo(prejudice,[2,1,0,0],'negative minus reversed positive; non-prejudiced'),
  prejudiced:linearCombo(prejudice,[2,1,2,1],'negative minus reversed positive; prejudiced'),
  difference:linearCombo(prejudice,[0,0,2,1],'prejudiced minus non-prejudiced asymmetry')};

function counts(data){const out={};for(const d of data){const k=d.highPrej?'prejudiced':'non_prejudiced';if(!out[k])out[k]={rows:0,negative:0,samples:new Set()};out[k].rows++;out[k].negative+=d.neg;out[k].samples.add(d.sample);}return Object.fromEntries(Object.entries(out).map(([k,v])=>[k,{rows:v.rows,negative:v.negative,negative_fraction:v.negative/v.rows,samples:v.samples.size}]));}
function samplePrejudiceBalance() {
  const m=new Map();
  for(const x of eligible){if(![1,2,3].includes(num(x.v8_Prej_Raw)))continue;const s=sampleKey(x);if(!m.has(s))m.set(s,{high:x.v8_Prej_Raw===1||x.v8_Prej_Raw===3,n:0,neg:0});const z=m.get(s);z.n++;z.neg+=String(x.v45_Recode_Pos_Ambiv_Neg_Unc)==='1'?1:0;}
  const o={};for(const z of m.values()){const k=z.high?'prejudiced':'non_prejudiced';if(!o[k])o[k]={samples:0,mean_negative_fraction:0};o[k].samples++;o[k].mean_negative_fraction+=z.neg/z.n;}for(const z of Object.values(o))z.mean_negative_fraction/=z.samples;return o;
}
const eligiblePrejudiceRows=eligible.filter(x=>[1,2,3].includes(num(x.v8_Prej_Raw))).map(x=>({highPrej:(x.v8_Prej_Raw===1||x.v8_Prej_Raw===3),neg:String(x.v45_Recode_Pos_Ambiv_Neg_Unc)==='1',sample:sampleKey(x)}));
const eligiblePrejudiceCounts={};for(const x of eligiblePrejudiceRows){const k=x.highPrej?'prejudiced':'non_prejudiced';if(!eligiblePrejudiceCounts[k])eligiblePrejudiceCounts[k]={rows:0,negative:0,samples:new Set()};eligiblePrejudiceCounts[k].rows++;eligiblePrejudiceCounts[k].negative+=x.neg?1:0;eligiblePrejudiceCounts[k].samples.add(x.sample);}for(const x of Object.values(eligiblePrejudiceCounts)){x.negative_fraction=x.negative/x.rows;x.samples=x.samples.size;}
const result={source_rows:rows.length,eligible_rows:eligible.length,matched_groups:matchedGroups.length,matched_samples:new Set(all.map(x=>x.sample)).size,all_rows:all.length,basic,basicAsymmetry,prejudice,asymmetry,prejudice_counts:counts(classified),eligible_prejudice_counts:eligiblePrejudiceCounts,sample_prejudice_balance:samplePrejudiceBalance(),notes:[
  'All-row model uses inverse Fisher-z sampling variances and sample-cluster-robust SEs; it is a practical multilevel/meta-regression sensitivity model, not a substitute for the original raw-data covariance structure.',
  'High prejudice is the workbook classification v8 raw = 1 or 3; raw = -99 is excluded from the prejudice analysis.',
  'The prejudice moderator is sample-level and categorical, not a continuous measure of prejudice strength.'
]};
fs.writeFileSync(`${root}/report/all_multilevel.json`,JSON.stringify(result,null,2));
console.log(JSON.stringify(result,null,2));
