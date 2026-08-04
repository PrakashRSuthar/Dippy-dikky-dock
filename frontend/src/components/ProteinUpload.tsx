// ProteinUpload.tsx
// When multiple ligands: asks user "Run together (parallel)" or "Run as batch (sequential)"
// Batch = one job per ligand queued one-by-one, each gets own job ID
// Together = all submitted at once via /api/dock/batch
import { useState } from 'react';
import { Settings, Activity, Upload, Atom, X, Layers, Zap } from 'lucide-react';
import { useApiHealth } from '../hooks/useApiHealth';
import { useFileUpload } from '../hooks/useFileUpload';
import { useDocking } from '../hooks/useDocking';
import { MoleculeVisualization } from './MoleculeVisualization';

interface ProteinUploadProps {
  onDockingStarted: (jobId: string, allJobIds?: string[], batchId?: string) => void;
}

const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

function saveBatchState(batchId: string, jobIds: string[], labels: Record<string, string>) {
  try {
    const raw = localStorage.getItem('dock_batches') || '{}';
    const map = JSON.parse(raw);
    map[batchId] = jobIds;
    localStorage.setItem('dock_batches', JSON.stringify(map));
    const existing = JSON.parse(localStorage.getItem('dock_job_labels') || '{}');
    localStorage.setItem('dock_job_labels', JSON.stringify({ ...existing, ...labels }));
  } catch {}
}

// ── Run mode dialog ───────────────────────────────────────────────────────────
interface RunModeDialogProps {
  count: number;
  onSelect: (mode: 'batch' | 'together') => void;
  onCancel: () => void;
}
function RunModeDialog({ count, onSelect, onCancel }: RunModeDialogProps) {
  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
      <div className="bg-white rounded-2xl shadow-2xl w-full max-w-md p-6">
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-lg font-semibold">Run {count} ligands — choose mode</h2>
          <button onClick={onCancel} className="text-gray-400 hover:text-gray-600"><X className="w-5 h-5"/></button>
        </div>
        <p className="text-sm text-gray-500 mb-6">
          You have <strong>{count} ligands</strong> selected against the same protein. How would you like to run them?
        </p>
        <div className="space-y-3">
          <button
            onClick={() => onSelect('batch')}
            className="w-full flex items-start gap-4 p-4 border-2 border-gray-200 rounded-xl hover:border-blue-500 hover:bg-blue-50 transition-colors text-left group"
          >
            <div className="w-10 h-10 bg-blue-100 rounded-lg flex items-center justify-center shrink-0 group-hover:bg-blue-200">
              <Layers className="w-5 h-5 text-blue-600"/>
            </div>
            <div>
              <div className="font-medium text-gray-900">Run as batch (one by one)</div>
              <div className="text-sm text-gray-500 mt-0.5">
                Each ligand gets its own job ID. Jobs run sequentially — safe for low-RAM systems. Results have separate tabs.
              </div>
            </div>
          </button>

          <button
            onClick={() => onSelect('together')}
            className="w-full flex items-start gap-4 p-4 border-2 border-gray-200 rounded-xl hover:border-green-500 hover:bg-green-50 transition-colors text-left group"
          >
            <div className="w-10 h-10 bg-green-100 rounded-lg flex items-center justify-center shrink-0 group-hover:bg-green-200">
              <Zap className="w-5 h-5 text-green-600"/>
            </div>
            <div>
              <div className="font-medium text-gray-900">Run together (parallel)</div>
              <div className="text-sm text-gray-500 mt-0.5">
                All ligands submitted at once. Faster on multi-core systems. Results appear in tabs as each completes.
              </div>
            </div>
          </button>
        </div>
        <button onClick={onCancel} className="w-full mt-4 py-2 text-sm text-gray-500 hover:text-gray-700">Cancel</button>
      </div>
    </div>
  );
}

// ── Main component ────────────────────────────────────────────────────────────
export const ProteinUpload = ({ onDockingStarted }: ProteinUploadProps) => {
  const [showSettings,      setShowSettings]      = useState(false);
  const [showRunModeDialog, setShowRunModeDialog] = useState(false);

  // Protein
  const [activeProteinTab, setActiveProteinTab] = useState<'database'|'upload'>('database');
  const [proteinId,        setProteinId]        = useState('1HSG');
  const [proteinInput,     setProteinInput]     = useState('');
  const [proteinSource,    setProteinSource]    = useState<'database'|'upload'|null>(null);

  // Ligands
  const [activeLigandTab, setActiveLigandTab] = useState<'database'|'upload'>('database');
  const [ligandId,        setLigandId]        = useState('aspirin');
  const [ligandInput,     setLigandInput]     = useState('');   // currently previewed
  const [_ligandSource,    _setLigandSource]    = useState<'database'|'upload'|null>(null);
  const [ligandInputs,    setLigandInputs]    = useState<string[]>([]);
  const [ligandFiles,     setLigandFiles]     = useState<{name:string;path:string}[]>([]);

  // Retention
  const [retentionPolicy, setRetentionPolicy] = useState<'save'|'delete'|'keep7d'>('keep7d');

  // Submitting state
  const [isSubmitting, setIsSubmitting] = useState(false);

  const { isHealthy, isChecking }                      = useApiHealth();
  const { uploadFile, isUploading }                    = useFileUpload();
  const { startDocking, isStarting }                   = useDocking();

  // ── Cleaning policy ──────────────────────────────────────────────────────
  const cleaningPolicy = { keep_waters:false, keep_ions:true, keep_solvents:false, keep_cofactors:true };

  // ── Protein handlers ─────────────────────────────────────────────────────
  const handleProteinUpload = async (file: File) => {
    if (!isHealthy) { alert('Backend offline. Start the FastAPI server first.'); return; }
    const result = await uploadFile(file, 'protein');
    if (result?.file_path) { setProteinInput(result.file_path); setProteinSource('upload'); }
  };

  const handleProteinDatabase = () => {
    if (proteinId.trim().length >= 4) { setProteinInput(proteinId.trim()); setProteinSource('database'); }
  };

  // ── Ligand handlers ──────────────────────────────────────────────────────
  const addLigandPaths = (paths: string[], names?: string[]) => {
    setLigandInputs(prev => {
      const seen = new Set(prev);
      const added = paths.filter(p => !seen.has(p));
      return [...prev, ...added];
    });
    setLigandFiles(prev => {
      const seen = new Set(prev.map(f => f.path));
      const added = paths.filter(p => !seen.has(p))
        .map((p, i) => ({ path: p, name: names?.[i] || p.split('/').pop() || p }));
      return [...prev, ...added];
    });
    setLigandInput(prev => prev || paths[0]);
    _setLigandSource('database');
  };

  const handleLigandDatabase = () => {
    const raw = ligandId.trim();
    if (!raw) return;
    const list = raw.split(',').map(s => s.trim()).filter(Boolean);
    if (!list.length) return;
    addLigandPaths(list, list);
  };

  const handleLigandFilesUpload = async (files: FileList | File[]) => {
    if (!isHealthy) { alert('Backend offline.'); return; }
    const accepted = Array.from(files).filter(f =>
      ['.sdf','.mol','.mol2','.pdbqt'].some(ext => f.name.toLowerCase().endsWith(ext)));
    if (!accepted.length) return;
    const uploaded: {name:string;path:string}[] = [];
    for (const file of accepted) {
      const result = await uploadFile(file, 'ligand');
      if (result?.file_path) uploaded.push({ name: file.name, path: result.file_path });
    }
    if (uploaded.length) {
      addLigandPaths(uploaded.map(u => u.path), uploaded.map(u => u.name));
      _setLigandSource('upload');
    }
  };

  const removeLigand = (idx: number) => {
    const removed = ligandInputs[idx];
    const nextInputs = ligandInputs.filter((_,i) => i !== idx);
    const nextFiles  = ligandFiles.filter((_,i) => i !== idx);
    setLigandInputs(nextInputs);
    setLigandFiles(nextFiles);
    if (ligandInput === removed) setLigandInput(nextInputs[0] || '');
  };

  // ── Docking submission ───────────────────────────────────────────────────
  const effectiveLigands = ligandInputs.length > 0 ? ligandInputs : (ligandInput ? [ligandInput] : []);

  const handleStartDockingClick = () => {
    if (!proteinInput || effectiveLigands.length === 0) return;
    if (effectiveLigands.length > 1) {
      setShowRunModeDialog(true);  // ask user
    } else {
      runSingle(effectiveLigands[0]);
    }
  };

  // Single job
  const runSingle = async (ligand: string) => {
    setIsSubmitting(true);
    const jobId = await startDocking({
      protein_input:   proteinInput,
      ligand_input:    ligand,
      job_name:        `${proteinInput}_${ligand}_${Date.now()}`,
      cleaning_policy: cleaningPolicy,
      retention:       retentionPolicy,
    });
    setIsSubmitting(false);
    if (jobId) onDockingStarted(jobId, [jobId]);
  };

  // Together: all submitted at once to /api/dock/batch
  const runTogether = async () => {
    setShowRunModeDialog(false);
    setIsSubmitting(true);
    try {
      const payload = {
        items: effectiveLigands.map((lig, i) => ({
          protein_input: proteinInput,
          ligand_input:  lig,
          job_name:      ligandFiles[i]?.name || lig,
        })),
        cleaning_policy: cleaningPolicy,
        retention:       retentionPolicy,
      };
      const res  = await fetch(`${apiBase}/api/dock/batch`, {
        method: 'POST', headers: {'Content-Type':'application/json'},
        body: JSON.stringify(payload),
      });
      if (!res.ok) throw new Error(await res.text());
      const data   = await res.json();
      const jobIds: string[] = (data.jobs || []).map((j: any) => j.job_id);
      const labels: Record<string,string> = {};
      (data.jobs || []).forEach((j: any, i: number) => {
        labels[j.job_id] = ligandFiles[i]?.name || effectiveLigands[i] || j.job_id.slice(0,8);
      });
      saveBatchState(data.batch_id, jobIds, labels);
      onDockingStarted(jobIds[0], jobIds, data.batch_id);
    } catch (e) {
      console.error(e); alert('Failed to start parallel docking: ' + e);
    }
    setIsSubmitting(false);
  };

  // Batch: one by one sequentially — submit all, track all job IDs
  const runBatch = async () => {
    setShowRunModeDialog(false);
    setIsSubmitting(true);
    const batchId = crypto.randomUUID();
    const jobIds:  string[] = [];
    const labels:  Record<string,string> = {};
    try {
      for (let i = 0; i < effectiveLigands.length; i++) {
        const lig   = effectiveLigands[i];
        const label = ligandFiles[i]?.name || lig;
        const res   = await fetch(`${apiBase}/api/dock`, {
          method: 'POST', headers: {'Content-Type':'application/json'},
          body: JSON.stringify({
            protein_input:   proteinInput,
            ligand_input:    lig,
            job_name:        label,
            cleaning_policy: cleaningPolicy,
            retention:       retentionPolicy,
          }),
        });
        if (!res.ok) { console.error(`Job ${i} failed to start`); continue; }
        const data = await res.json();
        if (data.job_id) { jobIds.push(data.job_id); labels[data.job_id] = label; }
      }
      if (!jobIds.length) { alert('No jobs started.'); setIsSubmitting(false); return; }
      saveBatchState(batchId, jobIds, labels);
      onDockingStarted(jobIds[0], jobIds, batchId);
    } catch (e) {
      console.error(e); alert('Batch submission error: ' + e);
    }
    setIsSubmitting(false);
  };

  const canProceed    = Boolean(proteinInput) && effectiveLigands.length > 0 && isHealthy;
  const isValidProtId = proteinId.trim().length >= 4;
  const busy          = isSubmitting || isStarting || isUploading;

  // ── Render ────────────────────────────────────────────────────────────────
  return (
    <div className="max-w-full mx-auto p-4 lg:p-6">

      {/* Run mode dialog */}
      {showRunModeDialog && (
        <RunModeDialog
          count={effectiveLigands.length}
          onSelect={mode => mode === 'batch' ? runBatch() : runTogether()}
          onCancel={() => setShowRunModeDialog(false)}
        />
      )}

      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div className="flex items-center gap-3">
          <div className="w-10 h-10 bg-blue-600 rounded-lg flex items-center justify-center">
            <Activity className="w-6 h-6 text-white"/>
          </div>
          <div>
            <h1 className="text-2xl font-bold">MolecularDock Pro</h1>
            <p className="text-sm text-gray-600">Professional Molecular Docking Platform</p>
          </div>
          <div className="flex items-center gap-2 ml-4">
            <div className={`w-3 h-3 rounded-full ${isChecking?'bg-yellow-400 animate-pulse':isHealthy?'bg-green-400':'bg-red-400'}`}/>
            <span className="text-xs text-gray-600">{isChecking?'Checking…':isHealthy?'System OK':'Backend Offline'}</span>
          </div>
        </div>
        <button onClick={() => setShowSettings(true)} className="flex items-center gap-2 px-4 py-2 text-gray-600 hover:bg-gray-100 rounded-lg">
          <Settings className="w-4 h-4"/><span>Settings</span>
        </button>
      </div>

      {!isHealthy && (
        <div className="mb-6 p-4 bg-yellow-50 border border-yellow-200 rounded-lg">
          <p className="text-yellow-800 font-medium text-sm">Backend required for docking</p>
          <p className="text-yellow-700 text-sm mt-1">
            Run: <code className="bg-yellow-100 px-1 rounded">python api/docking_api.py</code>
          </p>
        </div>
      )}

      <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">

        {/* ── Protein column ─────────────────────────────────────────────── */}
        <div className="space-y-6">
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center gap-2">
                🧬 Protein Structure {proteinInput && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-500 mt-0.5">Upload or fetch from database</p>
            </div>

            <div className="flex border-b">
              {(['database','upload'] as const).map(tab => (
                <button key={tab} onClick={() => setActiveProteinTab(tab)}
                  className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 capitalize ${
                    activeProteinTab===tab ? 'text-blue-600 border-blue-600 bg-blue-50' : 'text-gray-600 border-transparent'}`}>
                  {tab}
                </button>
              ))}
            </div>

            <div className="p-4">
              {activeProteinTab === 'database' ? (
                <div className="space-y-3">
                  <input type="text" value={proteinId} onChange={e => setProteinId(e.target.value.toUpperCase())}
                    placeholder="1HSG or P12345" onKeyDown={e => e.key==='Enter' && handleProteinDatabase()}
                    className="w-full px-3 py-2 border rounded-md text-center font-mono focus:ring-2 focus:ring-blue-500 focus:outline-none"/>
                  <p className="text-xs text-gray-500 text-center">PDB ID or UniProt ID</p>
                  <button onClick={handleProteinDatabase} disabled={!isValidProtId}
                    className={`w-full py-2 rounded-md font-medium ${isValidProtId?'bg-blue-600 text-white hover:bg-blue-700':'bg-gray-200 text-gray-400 cursor-not-allowed'}`}>
                    Use Protein ID
                  </button>
                </div>
              ) : (
                <div onDrop={e => { e.preventDefault(); const f=e.dataTransfer.files[0]; if(f) handleProteinUpload(f); }}
                  onDragOver={e => e.preventDefault()}
                  className="border-2 border-dashed rounded-xl p-6 text-center hover:border-blue-400 transition-colors">
                  <Upload className="w-6 h-6 mx-auto mb-2 text-gray-400"/>
                  <p className="text-sm text-gray-500 mb-2">Drop PDB / PDBQT here</p>
                  <input type="file" accept=".pdb,.pdbqt" className="hidden" id="prot-file"
                    onChange={e => { const f=e.target.files?.[0]; if(f) handleProteinUpload(f); }}/>
                  <label htmlFor="prot-file" className="cursor-pointer text-xs border border-gray-300 rounded px-3 py-1 hover:bg-gray-50">Browse</label>
                </div>
              )}
              {proteinInput && (
                <div className="mt-3 p-2 bg-green-50 rounded-lg text-sm text-green-800 flex items-center gap-2">
                  <span>✅</span>
                  <span className="truncate">{proteinSource==='upload' ? proteinInput.split('/').pop() : proteinInput}</span>
                </div>
              )}
            </div>
          </div>

          {/* Protein preview */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h3 className="text-lg font-semibold">🔬 Protein Preview</h3>
            </div>
            <div className="p-4">
              <div className="h-80 border border-gray-200 rounded-lg bg-gray-50">
                <MoleculeVisualization moleculePath={proteinInput||null} moleculeType="protein" height="100%"/>
              </div>
            </div>
          </div>
        </div>

        {/* ── Ligand column ──────────────────────────────────────────────── */}
        <div className="space-y-6">
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center gap-2">
                💊 Ligand{effectiveLigands.length > 1 ? 's' : ''}
                {effectiveLigands.length > 1 && (
                  <span className="text-xs bg-blue-100 text-blue-800 px-2 py-0.5 rounded-full font-normal">
                    {effectiveLigands.length} selected
                  </span>
                )}
                {effectiveLigands.length > 0 && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-500 mt-0.5">Upload files or fetch from PubChem</p>
            </div>

            <div className="flex border-b">
              {(['database','upload'] as const).map(tab => (
                <button key={tab} onClick={() => setActiveLigandTab(tab)}
                  className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 capitalize ${
                    activeLigandTab===tab ? 'text-green-600 border-green-600 bg-green-50' : 'text-gray-600 border-transparent'}`}>
                  {tab==='database' ? 'PubChem' : 'Upload'}
                </button>
              ))}
            </div>

            <div className="p-4">
              {activeLigandTab === 'database' ? (
                <div className="space-y-3">
                  <input type="text" value={ligandId} onChange={e => setLigandId(e.target.value)}
                    onKeyDown={e => e.key==='Enter' && handleLigandDatabase()}
                    placeholder="aspirin, ibuprofen, 2244 — comma-separated for multiple"
                    className="w-full px-3 py-2 border rounded-md text-center font-mono focus:ring-2 focus:ring-green-500 focus:outline-none text-sm"/>
                  <p className="text-xs text-gray-500 text-center">Compound name or PubChem CID · use commas for multiple</p>
                  <button onClick={handleLigandDatabase} disabled={ligandId.trim().length < 2}
                    className={`w-full py-2 rounded-md font-medium ${ligandId.trim().length>=2?'bg-green-600 text-white hover:bg-green-700':'bg-gray-200 text-gray-400 cursor-not-allowed'}`}>
                    Add Ligand{ligandId.includes(',') ? 's' : ''}
                  </button>
                </div>
              ) : (
                <div onDrop={e => { e.preventDefault(); handleLigandFilesUpload(e.dataTransfer.files); }}
                  onDragOver={e => e.preventDefault()}
                  className="border-2 border-dashed rounded-xl p-6 text-center hover:border-green-400 transition-colors">
                  <Atom className="w-6 h-6 mx-auto mb-2 text-gray-400"/>
                  <p className="text-sm text-gray-500 mb-2">Drop SDF / MOL / MOL2 / PDBQT here</p>
                  <input type="file" accept=".sdf,.mol,.mol2,.pdbqt" multiple className="hidden" id="lig-file"
                    onChange={e => { if(e.target.files) handleLigandFilesUpload(e.target.files); }}/>
                  <label htmlFor="lig-file" className="cursor-pointer text-xs border border-gray-300 rounded px-3 py-1 hover:bg-gray-50">Browse</label>
                  {isUploading && <p className="mt-2 text-xs text-blue-600">Uploading…</p>}
                </div>
              )}

              {/* Ligand list with remove */}
              {ligandFiles.length > 0 && (
                <div className="mt-4 space-y-1">
                  <p className="text-xs font-medium text-gray-600 mb-2">{ligandFiles.length} ligand{ligandFiles.length>1?'s':''} queued:</p>
                  <div className="max-h-40 overflow-y-auto space-y-1">
                    {ligandFiles.map((f, idx) => (
                      <div key={`${f.path}-${idx}`}
                        onClick={() => setLigandInput(f.path)}
                        className={`flex items-center justify-between px-3 py-2 rounded-lg border cursor-pointer text-sm transition-colors ${
                          ligandInput===f.path ? 'bg-green-50 border-green-300 text-green-800' : 'bg-gray-50 border-gray-200 hover:bg-gray-100'}`}>
                        <span className="truncate mr-2">{f.name}</span>
                        <button onClick={e => { e.stopPropagation(); removeLigand(idx); }}
                          className="text-gray-400 hover:text-red-500 shrink-0">
                          <X className="w-4 h-4"/>
                        </button>
                      </div>
                    ))}
                  </div>
                  <p className="text-xs text-gray-400 mt-1">Click a ligand to preview it below</p>
                </div>
              )}
            </div>
          </div>

          {/* Ligand preview */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h3 className="text-lg font-semibold flex items-center gap-2">
                ⚗️ Ligand Preview
                {ligandInput && <span className="text-sm font-normal text-gray-500">— {ligandFiles.find(f=>f.path===ligandInput)?.name || ligandInput}</span>}
              </h3>
            </div>
            <div className="p-4">
              <div className="h-80 border border-gray-200 rounded-lg bg-gray-50">
                <MoleculeVisualization
                  moleculeType="ligand"
                  moleculePath={ligandInput || null}
                  height="100%"
                  ligandList={ligandFiles}
                  onLigandChange={path => setLigandInput(path)}
                />
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* ── Retention policy ──────────────────────────────────────────────── */}
      {canProceed && (
        <div className="mt-8 bg-white rounded-xl border shadow-sm p-6">
          <h3 className="text-lg font-semibold mb-4">💾 Data Retention</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            {[
              { value:'save',   icon:'💾', label:'Save permanently', desc:'Move to data/results/ after completion' },
              { value:'keep7d', icon:'⏳', label:'Keep temporarily',  desc:'Auto-delete after 7 days' },
              { value:'delete', icon:'🗑️', label:'Delete immediately',desc:'Remove files once results are shown' },
            ].map(opt => (
              <label key={opt.value} className={`flex items-start gap-3 p-4 border-2 rounded-xl cursor-pointer transition-colors ${
                retentionPolicy===opt.value ? 'border-blue-500 bg-blue-50' : 'border-gray-200 hover:bg-gray-50'}`}>
                <input type="radio" name="retention" value={opt.value} checked={retentionPolicy===opt.value}
                  onChange={e => setRetentionPolicy(e.target.value as any)} className="mt-1"/>
                <div>
                  <div className="font-medium flex items-center gap-2"><span>{opt.icon}</span>{opt.label}</div>
                  <p className="text-sm text-gray-500 mt-0.5">{opt.desc}</p>
                </div>
              </label>
            ))}
          </div>
        </div>
      )}

      {/* ── Start button ──────────────────────────────────────────────────── */}
      {canProceed && (
        <div className="mt-8 flex flex-col items-center gap-3">
          {effectiveLigands.length > 1 && (
            <p className="text-sm text-gray-500">
              {effectiveLigands.length} ligands selected · you'll choose batch or parallel after clicking
            </p>
          )}
          <button onClick={handleStartDockingClick} disabled={!canProceed || busy}
            className={`px-10 py-4 rounded-xl font-semibold text-lg transition-colors flex items-center gap-3 shadow-lg ${
              canProceed && !busy ? 'bg-blue-600 text-white hover:bg-blue-700' : 'bg-gray-300 text-gray-500 cursor-not-allowed'}`}>
            {busy ? (
              <>
                <svg className="w-5 h-5 animate-spin" fill="none" viewBox="0 0 24 24">
                  <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"/>
                  <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v8z"/>
                </svg>
                Starting…
              </>
            ) : (
              <><span>⚡</span>Start Molecular Docking</>
            )}
          </button>
        </div>
      )}

      {/* ── Settings modal ─────────────────────────────────────────────────── */}
      {showSettings && (
        <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50">
          <div className="bg-white rounded-xl p-6 w-96 shadow-2xl">
            <h3 className="text-lg font-semibold mb-4">Settings</h3>
            <label className="block text-sm font-medium mb-1">API Base URL</label>
            <input type="url" defaultValue={apiBase} className="w-full px-3 py-2 border rounded-md mb-4 focus:outline-none focus:ring-2 focus:ring-blue-500"/>
            <div className="flex gap-3">
              <button onClick={() => setShowSettings(false)} className="flex-1 py-2 border rounded-md hover:bg-gray-50">Cancel</button>
              <button onClick={() => setShowSettings(false)} className="flex-1 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700">Save</button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};
export default ProteinUpload;