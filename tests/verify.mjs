// No dependencies: execute with node tests/verify.mjs.
import fs from 'node:fs';
import vm from 'node:vm';
import assert from 'node:assert/strict';
const html=fs.readFileSync(new URL('../index.html',import.meta.url),'utf8');
const math=html.slice(html.indexOf('function qwp('),html.indexOf('function relativeOf('));
const context=vm.createContext({Math,Error,Number,Infinity});
vm.runInContext('const RAD=Math.PI/180;'+math+'; this.api={qwp,stokes,feature,leastSquares,norm2,condition};',context);
const {qwp,stokes,feature,leastSquares,norm2}=context.api;
const close=(a,b,t=1e-10)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
let tests=0;
for(const [theta,h,v] of [[0,1,0],[22.5,.75,.25],[45,.5,.5],[67.5,.75,.25],[90,1,0]]){
 const j=qwp(theta);close(norm2(j[0]),h);close(norm2(j[1]),v);tests++;
}
for(let angle=0;angle<180;angle+=.5){const a=qwp(angle),b=qwp(angle+180);a.flat().forEach((v,i)=>close(v,b.flat()[i]));close(norm2(a[0])+norm2(a[1]),1);tests++;}
// Independent Jones-generated observations validate the Mueller analyzer signs.
for(let i=0;i<100;i++){
 const amp=(i+.3)/101,delta=i*.713,j=[[Math.sqrt(amp),0],[Math.sqrt(1-amp)*Math.cos(delta),Math.sqrt(1-amp)*Math.sin(delta)]],s=stokes(j);
 for(const theta of [0,22.5,45,67.5,89,135,281]){const o=qwp(theta,j),d=(norm2(o[0])-norm2(o[1]))/(norm2(o[0])+norm2(o[1]));close(d,feature(theta).reduce((v,x,k)=>v+x*s[k],0));tests++;}
 const angles=[0,22.5,45,67.5],ds=angles.map(t=>{const o=qwp(t,j);return norm2(o[0])-norm2(o[1]);});
 [ds[0],ds[1]-ds[3],ds[2]].forEach((x,k)=>close(x,s[k]));
 const A=[0,22.5,45,67.5,91,143,270].map(t=>feature(t));const fit=leastSquares(A,A.map(a=>a.reduce((v,x,k)=>v+x*s[k],0)));fit.s.forEach((x,k)=>close(x,s[k]));tests++;
}
for(const s of [[1,0,0],[0,1,0],[0,0,1],[0,0,-1],[.36,.48,.8],[.1,-.2,.3],[0,0,0]]){const A=[0,22.5,45,67.5,131].map(t=>feature(t));const fit=leastSquares(A,A.map(a=>a.reduce((v,x,k)=>v+x*s[k],0)));fit.s.forEach((v,i)=>close(v,s[i]));tests++;}
for(const angles of [[0,90,180],[45,45,45,45],[0,0,0],[0,22.5]]){assert.throws(()=>leastSquares(angles.map(a=>feature(a)),angles.map(()=>0)));tests++;}
const A=[0,.00000001,.00000002].map(t=>feature(t));assert.throws(()=>leastSquares(A,[1,1,1]));tests++;
close(stokes(qwp(45))[2],-1);close(stokes(qwp(135))[2],1);tests+=2;
console.log(`PASS ${tests} physics checks: Jones/PBS, 180° periodicity, Stokes signs, basic inversion, QR, rank/conditioning.`);
// Data handling checks reuse the production parser and estimator in an isolated VM.
const helpers=html.slice(html.indexOf('function number('),html.indexOf('function fmt('));
const validation=html.slice(html.indexOf('function validateState('),html.indexOf('try{const saved='));
const estimator=html.slice(html.indexOf('function relativeOf('),html.indexOf('function powerMarkup('));
const csv=html.slice(html.indexOf('const CSV_HEADER='),html.indexOf("$('export-csv').onclick="));
vm.runInContext(`const BASE=[0,22.5,45,67.5];const wrapAngle=n=>((n%360)+360)%360;let state;function trimNum(n,d=3){return Number.isFinite(n)?String(Number(n.toFixed(d))):'—';}${helpers}${validation}${estimator}${csv};this.data={set:s=>state=s,get:()=>state,estimate,csvText,importCSV,parseCSV,number};`,context);
const data=context.data;
const original={version:1,unit:'mW',theta0:'350',mapping:'AV',etaA:'.8',etaB:'.6',convention:'minus',method:'ls',rows:[]};
const target=[.36,.48,-.8];
for(const [i,theta] of [0,22.5,45,67.5,131].entries()){const d=feature(theta,-1).reduce((v,x,k)=>v+x*target[k],0);original.rows.push({kind:i<4?'base':'extra',...(i<4?{relative:theta}:{dial:'121'}),a:String(.05*(1-d)*.8),b:String(.05*(1+d)*.6)});}
data.set(original);let estimate=data.estimate();estimate.s.forEach((v,i)=>close(v,target[i]));close(estimate.S0,.1);
const exported=data.csvText(),imported=data.importCSV(exported);assert.deepEqual(JSON.parse(JSON.stringify(imported)),original);data.set(imported);estimate=data.estimate();estimate.s.forEach((v,i)=>close(v,target[i]));
assert.throws(()=>data.importCSV(exported.replace(',AV,',',INVALID,')));assert.throws(()=>data.importCSV(exported.replace('0.0256','-2')));assert.throws(()=>data.parseCSV('"unterminated'));
assert.deepEqual(JSON.parse(JSON.stringify(data.parseCSV('a,b\r\n"hello, world","quote ""x"""'))),[['a','b'],['hello, world','quote "x"']]);
for(const x of ['',' ','NaN','Infinity','1;alert(1)','0xFF','1,2'])assert.ok(Number.isNaN(data.number(x)));
for(const [s,v] of [['.5',.5],['1e-3',.001],['-3',-3]])close(data.number(s),v);
original.method='basic';original.unit='uW';original.theta0='0';original.mapping='AH';original.etaA='1';original.etaB='1';original.convention='plus';original.rows=[100,100,100,0].map((h,i)=>({kind:'base',relative:[0,22.5,45,67.5][i],a:String(h),b:String(100-h)}));data.set(original);assert.ok(data.estimate().p>1);original.rows[0].a='-1';assert.throws(()=>data.estimate());original.rows[0].a='0';original.rows[0].b='0';assert.throws(()=>data.estimate());
console.log('PASS data checks: CSV round trip, malformed CSV, calibrated estimator, numeric validation, invalid power, DOP not clamped.');
