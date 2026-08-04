// ResultsPage.tsx
// Per-ligand tabs: each tab = one job ID, shows live status badge, clicking loads that job's results
import { useState, useEffect, useMemo, useCallback } from 'react';
import { ArrowLeft, Download, RefreshCw, BarChart3, Zap, Atom, CheckCircle,
         XCircle, FileIcon, CheckCircle2, Circle, XCircle as XIcon } from 'lucide-react';
import { MoleculeVisualization } from '../components/MoleculeVisualization';
import { ResultsTable }          from '../components/ResultsTable';
import { ResultsChart }          from '../components/ResultsChart';

interface ResultsPageProps {
  jobId?:   string;
  jobIds?:  string[];
  batchId?: string;
  onBack:   () => void;
}

interface PocketInfo { center: number[]; size: number[]; confidence: string; method: string }
interface FileInfo   { name: string; filename: string; path: string; size: string; download_url: string }
interface PoseRow    { pose: number; affinity: string; rmsd_lb: string; rmsd_ub: string }
interface DockingResults {
  job_id: string; status: string; best_affinity: number; total_poses: number;
  pocket_center: number[]; pocket_size: number[]; method: string; confidence: string;
  run_directory: string; all_pockets: PocketInfo[]; files: Record<string,string>;
  protein_analysis: string; cleaning_policy: any; retention: string;
  drift_events?: number; consensus_clusters?: number; methods_used?: string[];
  pose_data?: { pose:number; affinity:number; rmsd_lb:number; rmsd_ub:number }[];
}
type JobStatus = 'queued'|'running'|'completed'|'failed'|'unknown';

const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

function loadBatchJobs(batchId?: string): string[] | null {
  if (!batchId) return null;
  try {
    const raw = localStorage.getItem('dock_batches');
    if (!raw) return null;
    const map = JSON.parse(raw) as Record<string,string[]>;
    return Array.isArray(map[batchId]) && map[batchId].length ? map[batchId] : null;
  } catch { return null; }
}
function loadJobLabels(): Record<string,string> {
  try { return JSON.parse(localStorage.getItem('dock_job_labels') || '{}'); }
  catch { return {}; }
}

// ── Status badge ──────────────────────────────────────────────────────────────
function StatusBadge({ st }: { st: JobStatus }) {
  const base = 'inline-flex items-center gap-1 px-2 py-0.5 rounded-full text-xs font-medium';
  if (st==='completed') return <span className={`${base} bg-green-100 text-green-700`}><CheckCircle2 className="w-3 h-3"/>Done</span>;
  if (st==='running')   return <span className={`${base} bg-blue-100 text-blue-700`}><Circle className="w-3 h-3 animate-pulse"/>Running</span>;
  if (st==='failed')    return <span className={`${base} bg-red-100 text-red-700`}><XIcon className="w-3 h-3"/>Failed</span>;
  if (st==='queued')    return <span className={`${base} bg-yellow-100 text-yellow-700`}><Circle className="w-3 h-3"/>Queued</span>;
  return <span className={`${base} bg-gray-100 text-gray-600`}>—</span>;
}

export const ResultsPage = ({ jobId, jobIds, batchId, onBack }: ResultsPageProps) => {
  const labels = useMemo(() => loadJobLabels(), []);

  // Build job list from all possible sources
  const allJobIds = useMemo<string[]>(() => {
    if (jobIds?.length) return jobIds;
    const rec = loadBatchJobs(batchId);
    if (rec?.length) return rec;
    return jobId ? [jobId] : [];
  }, [jobId, jobIds, batchId]);

  const [jobsList,     setJobsList]     = useState(allJobIds);
  const [activeIdx,    setActiveIdx]    = useState(0);
  const [jobStatuses,  setJobStatuses]  = useState<Record<string,JobStatus>>({});
  const [results,      setResults]      = useState<DockingResults|null>(null);
  const [files,        setFiles]        = useState<FileInfo[]>([]);
  const [poseData,     setPoseData]     = useState<PoseRow[]>([]);
  const [selectedPose, setSelectedPose] = useState(0);
  const [loading,      setLoading]      = useState(true);
  const [error,        setError]        = useState<string|null>(null);
  const [resultsCache, setResultsCache] = useState<Record<string, DockingResults>>({});

  const activeJobId = jobsList[activeIdx] || '';

  // Keep list in sync
  useEffect(() => { if (allJobIds.length) setJobsList(allJobIds); }, [allJobIds]);

  // ── Poll all job statuses ────────────────────────────────────────────────
  useEffect(() => {
    if (!jobsList.length) return;
    let cancelled = false;
    const tick = async () => {
      const updates: Record<string,JobStatus> = {};
      await Promise.all(jobsList.map(async jid => {
        try {
          const r = await fetch(`${apiBase}/api/jobs/${jid}/status`);
          if (!r.ok) { updates[jid]='unknown'; return; }
          const j = await r.json();
          updates[jid] = (j?.status as JobStatus) || 'unknown';
          // Auto-load result into cache when newly completed
          if (j?.status === 'completed' && !resultsCache[jid]) {
            fetchAndCacheResult(jid);
          }
        } catch { updates[jid] = 'unknown'; }
      }));
      if (!cancelled) setJobStatuses(p => ({...p,...updates}));
      if (!cancelled) setTimeout(tick, 1500);
    };
    tick();
    return () => { cancelled = true; };
  }, [jobsList]);

  // ── Fetch + cache one job's results ──────────────────────────────────────
  const fetchAndCacheResult = useCallback(async (jid: string) => {
    try {
      const res = await fetch(`${apiBase}/api/jobs/${jid}/result`);
      if (!res.ok) return;
      const data: DockingResults = await res.json();
      setResultsCache(c => ({ ...c, [jid]: data }));
    } catch {}
  }, []);

  // ── Load active job's results ────────────────────────────────────────────
  const loadJob = useCallback(async (jid: string) => {
    if (!jid) return;
    setLoading(true); setError(null); setResults(null); setPoseData([]); setFiles([]);

    // Use cache if available
    if (resultsCache[jid]) {
      applyResult(resultsCache[jid]); setLoading(false); return;
    }

    try {
      const res = await fetch(`${apiBase}/api/jobs/${jid}/result`);
      if (!res.ok) {
        const err = await res.json().catch(() => ({}));
        throw new Error(err.detail || `HTTP ${res.status}`);
      }
      const data: DockingResults = await res.json();
      setResultsCache(c => ({ ...c, [jid]: data }));
      applyResult(data);
    } catch (e: any) {
      setError(e?.message || 'Unknown error');
    } finally {
      setLoading(false);
    }

    // Load files list
    try {
      const fr = await fetch(`${apiBase}/api/jobs/${jid}/files`);
      if (fr.ok) setFiles((await fr.json()).files || []);
    } catch {}
  }, [resultsCache]);

  const applyResult = (data: DockingResults) => {
    setResults(data);
    if (data.pose_data?.length) {
      setPoseData(data.pose_data.map(p => ({
        pose:     p.pose,
        affinity: p.affinity.toFixed(2),
        rmsd_lb:  p.rmsd_lb.toFixed(2),
        rmsd_ub:  p.rmsd_ub.toFixed(2),
      })));
    } else {
      setPoseData([{ pose:1, affinity: data.best_affinity.toFixed(2), rmsd_lb:'0.00', rmsd_ub:'0.00' }]);
    }
    setSelectedPose(0);
  };

  // Load when active tab changes
  useEffect(() => { if (activeJobId) loadJob(activeJobId); }, [activeJobId]);

  const downloadFile = (url: string, name: string) => {
    const a = document.createElement('a');
    a.href=`${apiBase}${url}`; a.download=name;
    document.body.appendChild(a); a.click(); document.body.removeChild(a);
  };

  const getTabLabel = (jid: string, idx: number) =>
    labels[jid] || `Ligand ${idx+1} (${jid.slice(0,6)}…)`;

  // ── States ────────────────────────────────────────────────────────────────
  const isBatch = jobsList.length > 1;

  if (!activeJobId) return (
    <div className="min-h-screen flex items-center justify-center">
      <p className="text-gray-500">No job selected.</p>
    </div>
  );

  return (
    <div className="min-h-screen bg-gray-50">

      {/* ── Header ─────────────────────────────────────────────────────────── */}
      <div className="bg-white border-b sticky top-0 z-10 shadow-sm">
        <div className="max-w-full mx-auto px-6 py-3">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-4">
              <button onClick={onBack} className="flex items-center gap-2 text-gray-600 hover:text-gray-900">
                <ArrowLeft className="w-5 h-5"/><span>Back</span>
              </button>
              <div>
                <h1 className="text-xl font-bold flex items-center gap-2">
                  <CheckCircle className="w-5 h-5 text-green-600"/>
                  {isBatch ? `Batch Results (${jobsList.length} ligands)` : 'Docking Results'}
                </h1>
                {results && <p className="text-xs text-gray-500 font-mono mt-0.5">Job: {results.job_id}</p>}
              </div>
            </div>
            {files.length > 0 && (
              <button onClick={() => files.forEach((f,i) => setTimeout(()=>downloadFile(f.download_url,f.filename||f.name),i*200))}
                className="flex items-center gap-2 px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700 text-sm">
                <Download className="w-4 h-4"/>Download All
              </button>
            )}
          </div>

          {/* ── Per-ligand tabs ─────────────────────────────────────────────── */}
          {isBatch && (
            <div className="mt-3 flex items-center gap-1 overflow-x-auto pb-1">
              {jobsList.map((jid, idx) => {
                const active = idx === activeIdx;
                const st     = jobStatuses[jid] || 'unknown';
                return (
                  <button key={jid} onClick={() => setActiveIdx(idx)} title={jid}
                    className={`shrink-0 flex items-center gap-2 px-3 py-2 rounded-lg border text-sm transition-all ${
                      active
                        ? 'bg-white border-blue-500 text-blue-700 shadow font-medium'
                        : 'bg-gray-50 border-gray-200 text-gray-600 hover:bg-white hover:border-gray-300'
                    }`}>
                    <span className="max-w-[120px] truncate">{getTabLabel(jid, idx)}</span>
                    <StatusBadge st={st}/>
                  </button>
                );
              })}
            </div>
          )}
        </div>
      </div>

      {/* ── Content ────────────────────────────────────────────────────────── */}
      <div className="max-w-full mx-auto px-4 py-6">

        {/* Loading state */}
        {loading && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 animate-spin text-blue-600 mx-auto mb-4"/>
              <p className="text-gray-600 text-sm">
                {jobStatuses[activeJobId]==='running' || jobStatuses[activeJobId]==='queued'
                  ? 'Job still running — waiting for results…'
                  : 'Loading results…'}
              </p>
              {(jobStatuses[activeJobId]==='running'||jobStatuses[activeJobId]==='queued') && (
                <p className="text-xs text-gray-400 mt-2">Results will appear automatically when complete.</p>
              )}
            </div>
          </div>
        )}

        {/* Error state */}
        {!loading && error && (
          <div className="flex items-center justify-center py-16">
            <div className="text-center bg-white p-8 rounded-xl shadow border border-red-200 max-w-md">
              <XCircle className="w-12 h-12 text-red-500 mx-auto mb-4"/>
              <h2 className="text-lg font-bold mb-2">Error Loading Results</h2>
              <p className="text-gray-600 text-sm mb-6">{error}</p>
              <div className="flex gap-3 justify-center">
                <button onClick={onBack} className="px-4 py-2 bg-gray-600 text-white rounded-lg hover:bg-gray-700">Back</button>
                <button onClick={() => loadJob(activeJobId)} className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700">Retry</button>
              </div>
            </div>
          </div>
        )}

        {/* Results */}
        {!loading && !error && results && (
          <div className="space-y-6">
            {/* Summary cards */}
            <div className="grid grid-cols-2 md:grid-cols-3 xl:grid-cols-5 gap-4">
              {[
                { label:'Best Affinity', val:`${results.best_affinity} kcal/mol`, color:'text-green-600', icon:<Zap className="w-7 h-7 text-green-500 opacity-60"/> },
                { label:'Poses Found',   val:results.total_poses,   color:'text-blue-600',   icon:<Atom className="w-7 h-7 text-blue-500 opacity-60"/> },
                { label:'Binding Sites', val:results.all_pockets?.length||1, color:'text-purple-600', icon:<BarChart3 className="w-7 h-7 text-purple-500 opacity-60"/> },
                { label:'Box Drift Fixed',val:results.drift_events??0, color:'text-orange-600', icon:<span className="text-2xl">🎯</span> },
                { label:'Methods Used',  val:(results.methods_used||[]).length, color:'text-teal-600',  icon:<span className="text-2xl">🔬</span> },
              ].map(({ label, val, color, icon }) => (
                <div key={label} className="bg-white rounded-xl border shadow-sm p-4">
                  <div className="flex items-center justify-between">
                    <div>
                      <p className="text-xs text-gray-500 mb-1">{label}</p>
                      <p className={`text-xl font-bold ${color} break-all leading-tight`}>{String(val)}</p>
                    </div>
                    {icon}
                  </div>
                </div>
              ))}
            </div>

            {/* Pocket methods used */}
            {results.methods_used?.length && (
              <div className="flex items-center gap-2 flex-wrap">
                <span className="text-sm text-gray-500">Pocket methods:</span>
                {results.methods_used.map(m => (
                  <span key={m} className="text-xs px-2 py-0.5 bg-blue-100 text-blue-800 rounded-full">{m}</span>
                ))}
                {results.consensus_clusters != null && (
                  <span className="text-xs text-gray-400">· {results.consensus_clusters} consensus clusters</span>
                )}
              </div>
            )}

            {/* 3D + Files/Pockets */}
            <div className="grid grid-cols-1 xl:grid-cols-3 gap-6">
              <div className="xl:col-span-2 bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b flex items-center gap-3">
                  <h2 className="text-lg font-semibold">🔬 Docked Complex</h2>
                  <span className="text-sm bg-green-100 text-green-800 px-2 py-0.5 rounded-full">Pose {selectedPose+1}</span>
                </div>
                <div className="p-4">
                  <div className="h-[50vh] min-h-[420px] border border-gray-200 rounded-lg bg-gray-50">
                    <MoleculeVisualization
                      proteinPath={results.files?.prepared_protein||results.files?.raw_protein}
                      dockedPath={results.files?.docked_poses}
                      height="100%" pocketInfo={{ center:results.pocket_center, size:results.pocket_size }}
                      showPockets isResultsMode selectedPose={selectedPose}
                      poseData={poseData.map(p=>({affinity:p.affinity}))}/>
                  </div>
                </div>
              </div>

              <div className="bg-white rounded-xl border shadow-sm flex flex-col overflow-hidden">
                {/* Files */}
                <div className="px-5 py-3 border-b"><h2 className="font-semibold">📁 Files</h2></div>
                <div className="p-4 flex-1 overflow-y-auto">
                  {files.length === 0
                    ? <p className="text-sm text-gray-400">No files available.</p>
                    : <ul className="space-y-2">
                        {files.map(f => (
                          <li key={f.path} className="flex items-center justify-between gap-2 p-2 border rounded-lg text-sm">
                            <div className="flex items-center gap-2 min-w-0">
                              <FileIcon className="w-4 h-4 text-gray-400 shrink-0"/>
                              <div className="min-w-0">
                                <p className="truncate font-medium text-xs">{f.filename||f.name}</p>
                                {f.size && <p className="text-xs text-gray-400">{f.size}</p>}
                              </div>
                            </div>
                            <button onClick={() => downloadFile(f.download_url, f.filename||f.name)}
                              className="shrink-0 flex items-center gap-1 px-2 py-1 text-xs bg-gray-100 hover:bg-gray-200 rounded">
                              <Download className="w-3 h-3"/>Get
                            </button>
                          </li>
                        ))}
                      </ul>}
                </div>

                {/* Pockets */}
                <div className="px-5 py-3 border-t border-b">
                  <h2 className="font-semibold flex items-center gap-2">
                    🎯 Binding Pockets
                    <span className="text-xs bg-purple-100 text-purple-800 px-2 py-0.5 rounded-full">{results.all_pockets?.length||0}</span>
                  </h2>
                </div>
                <div className="p-4 overflow-y-auto max-h-64 space-y-2">
                  {results.all_pockets?.map((p,i) => (
                    <div key={i} className={`p-3 border rounded-lg text-xs ${i===0?'bg-blue-50 border-blue-200':''}`}>
                      <div className="flex items-center justify-between mb-1">
                        <span className="font-medium text-sm">{i===0&&'★ '}Pocket {i+1}</span>
                        <span className={`px-1.5 py-0.5 rounded text-xs ${
                          p.confidence==='high'?'bg-green-100 text-green-800':
                          p.confidence==='medium'?'bg-yellow-100 text-yellow-800':'bg-gray-100 text-gray-600'}`}>
                          {p.confidence}
                        </span>
                      </div>
                      <p className="text-gray-600">Center: ({p.center.map(c=>c.toFixed(1)).join(', ')})</p>
                      <p className="text-gray-600">Size: {p.size.map(s=>s.toFixed(1)).join('×')} Å</p>
                      <p className="text-gray-500">Method: {p.method}</p>
                    </div>
                  ))}
                </div>
              </div>
            </div>

            {/* Chart + Table */}
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b"><h2 className="font-semibold">📊 Affinity Distribution</h2></div>
                <div className="p-6 h-64"><ResultsChart data={poseData} selectedPose={selectedPose}/></div>
              </div>
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="font-semibold">📋 Docking Poses
                    <span className="text-xs text-gray-400 font-normal ml-2">RMSD values from Vina output</span>
                  </h2>
                </div>
                <div className="p-4">
                  <ResultsTable data={poseData} selectedPose={selectedPose} onPoseSelect={setSelectedPose}/>
                </div>
              </div>
            </div>

            {/* Protein analysis */}
            {results.protein_analysis && (
              <div className="bg-white rounded-xl border shadow-sm p-4">
                <h3 className="font-semibold mb-2 text-sm">🧪 Protein Analysis</h3>
                <p className="text-sm text-gray-600">{results.protein_analysis}</p>
              </div>
            )}
          </div>
        )}

        {/* ── Queued / running state for batch ─────────────────────────────── */}
        {!loading && !error && !results && isBatch && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 animate-spin text-blue-600 mx-auto mb-4"/>
              <p className="text-gray-600">
                Job {activeIdx+1} is {jobStatuses[activeJobId]||'queued'} — select another tab or wait.
              </p>
            </div>
          </div>
        )}
      </div>
    </div>
  );
};
export default ResultsPage;