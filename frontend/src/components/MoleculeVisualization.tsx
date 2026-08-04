// MoleculeVisualization.tsx
// Fix: ligand preview now works for PubChem names, CIDs, file paths, SDF/PDBQT uploads
import React, { useEffect, useRef, useState, useCallback } from 'react';

declare global { interface Window { $3Dmol?: any } }

type MolFormat = 'pdb' | 'pdbqt' | 'sdf' | 'mol2';
interface Pocket  { center: number[]; size: number[] }
interface PoseRow { affinity: string }

interface MoleculeVisualizationProps {
  moleculePath?:  string | null;
  moleculeType?:  'protein' | 'ligand';
  proteinPath?:   string | null;
  dockedPath?:    string | null;
  height?:        string;
  pocketInfo?:    Pocket;
  showPockets?:   boolean;
  selectedPose?:  number;
  isResultsMode?: boolean;
  poseData?:      PoseRow[];
  ligandList?:    { name?: string; path: string }[];
  onLigandChange?: (path: string) => void;
  onProceed?:     (mode: 'current' | 'all', currentLigPath?: string) => void;
}

// ── Toolbar icon SVGs ────────────────────────────────────────────────────────
const Btn: React.FC<{ title: string; onClick: () => void; children: React.ReactNode }> = ({ title, onClick, children }) => (
  <button onClick={onClick} title={title}
    className="h-8 w-8 flex items-center justify-center rounded-md text-gray-600 hover:bg-gray-200 transition-colors">
    {children}
  </button>
);
const StyleSVG = () => (
  <svg className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
    <path d="M17.414 2.586a2 2 0 00-2.828 0L7 10.172V13h2.828l7.586-7.586a2 2 0 000-2.828z"/>
    <path fillRule="evenodd" d="M2 6a2 2 0 012-2h4a1 1 0 010 2H4v10h10v-4a1 1 0 112 0v4a2 2 0 01-2 2H4a2 2 0 01-2-2V6z" clipRule="evenodd"/>
  </svg>
);
const OptSVG = () => (
  <svg className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor">
    <path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.532 1.532 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.532 1.532 0 01-.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd"/>
  </svg>
);

export const MoleculeVisualization: React.FC<MoleculeVisualizationProps> = ({
  moleculePath, moleculeType = 'protein', proteinPath, dockedPath,
  height = '500px', pocketInfo, showPockets = true,
  isResultsMode = false, selectedPose = 0, poseData,
  ligandList, onLigandChange, onProceed,
}) => {
  const containerRef  = useRef<HTMLDivElement>(null);
  const viewerRef     = useRef<HTMLDivElement>(null);
  const viewerInst    = useRef<any>(null);
  const modelsLoaded  = useRef(false);

  const [isLoading,      setIsLoading]      = useState(false);
  const [error,          setError]          = useState<string | null>(null);
  const [dataSource,     setDataSource]     = useState('');
  const [openMenu,       setOpenMenu]       = useState<string | null>(null);
  const [surfaceVisible, setSurfaceVisible] = useState(false); // off by default — faster load
  const [labelsVisible,  setLabelsVisible]  = useState(false);
  const [pocketsVisible, setPocketsVisible] = useState(showPockets);
  const [isSpinning,     setIsSpinning]     = useState(false);
  const [ligandStyle,    setLigandStyle]    = useState<'ballstick'|'stick'|'line'>('ballstick');
  const [receptorStyle,  setReceptorStyle]  = useState<'cartoon'|'sticks'|'lines'>('cartoon');
  const [receptorColor,  setReceptorColor]  = useState<'lightblue'|'spectrum'|'white'|'grey'>('lightblue');
  const [ligandColor,    setLigandColor]    = useState<'Jmol'|'cpk'|'grey'>('Jmol');
  const [currentLig,     setCurrentLig]     = useState<string | undefined>(
    ligandList?.length ? ligandList[0].path : undefined
  );

  const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

  // ── Core: resolve any input string → { data, format, source } ───────────
  const resolveStructureData = useCallback(async (input: string) => {
    if (!input?.trim()) throw new Error('No input provided');
    const s = input.trim();

    // 1) Server file path or /api/download URL
    if (s.includes('/') || s.includes('\\')) {
      const url = s.startsWith('/api/download')
        ? `${apiBase}${s}`
        : `${apiBase}/api/download?path=${encodeURIComponent(s)}`;
      const res = await fetch(url);
      if (!res.ok) throw new Error(`Server file download failed (${res.status})`);
      const text = await res.text();
      if (!text || text.length < 10) throw new Error('Downloaded file is empty');
      const ext = (s.split('.').pop() || '').toLowerCase();
      const fmt = ['pdb','pdbqt','sdf','mol2'].includes(ext) ? ext as MolFormat : 'pdb';
      return { data: text, format: fmt, source: 'Server File' };
    }

    // 2) PDB ID (4 alphanumeric chars starting with digit)
    if (/^[1-9][A-Za-z0-9]{3}$/i.test(s)) {
      const res = await fetch(`https://files.rcsb.org/download/${s.toUpperCase()}.pdb`);
      if (res.ok) return { data: await res.text(), format: 'pdb', source: 'RCSB PDB' };
    }

    // 3) PubChem — try name then CID via multiple proxies
    const isNumeric = /^\d+$/.test(s);
    const searchType = isNumeric ? 'cid' : 'name';
    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${searchType}/${encodeURIComponent(s)}/SDF`;

    // Proxy 1: allorigins
    try {
      const res = await fetch(`https://api.allorigins.win/raw?url=${encodeURIComponent(pubchemUrl)}`);
      if (res.ok) {
        const text = await res.text();
        if (text && !text.includes('Status: 404') && text.includes('$$$$')) {
          return { data: text, format: 'sdf', source: `PubChem (${s})` };
        }
      }
    } catch { /* try next */ }

    // Proxy 2: corsproxy.io
    try {
      const res = await fetch(`https://corsproxy.io/?${encodeURIComponent(pubchemUrl)}`);
      if (res.ok) {
        const text = await res.text();
        if (text && !text.includes('Status: 404') && text.includes('$$$$')) {
          return { data: text, format: 'sdf', source: `PubChem (${s})` };
        }
      }
    } catch { /* give up */ }

    throw new Error(`Could not load structure for "${s}". Try a PDB ID or upload a file.`);
  }, [apiBase]);

  // ── Apply styles after model load ────────────────────────────────────────
  const applyStyles = useCallback(() => {
    if (!viewerInst.current || !modelsLoaded.current) return;
    const v = viewerInst.current;

    // Receptor
    const col = { color: receptorColor === 'spectrum' ? 'spectrum' : receptorColor };
    if (receptorStyle === 'cartoon')  v.setStyle({ hetflag: false }, { cartoon: { ...col, opacity: 0.9 } });
    else if (receptorStyle === 'sticks') v.setStyle({ hetflag: false }, { stick: { radius: 0.2, ...col } });
    else v.setStyle({ hetflag: false }, { line: col });

    // Ligand
    const scheme = ligandColor === 'grey' ? 'greyCarbon' : ligandColor;
    const tgt = isResultsMode ? { model: 1 } : (moleculeType === 'ligand' ? {} : {});
    if (ligandStyle === 'ballstick')
      v.setStyle(tgt, { stick: { radius: 0.3, colorscheme: scheme }, sphere: { radius: 0.5, colorscheme: scheme } });
    else if (ligandStyle === 'stick')
      v.setStyle(tgt, { stick: { radius: 0.25, colorscheme: scheme } });
    else
      v.setStyle(tgt, { line: { colorscheme: scheme } });

    // Surface (skip for ligand-only views — slow and meaningless)
    v.removeAllSurfaces();
    if (surfaceVisible && (moleculeType === 'protein' || isResultsMode)) {
      v.addSurface('MS', { opacity: 0.5, color: 'white' }, { model: 0 });
    }

    // Pocket box
    v.removeAllShapes(); v.removeAllLabels();
    if (pocketsVisible && pocketInfo) {
      const [x,y,z]    = pocketInfo.center;
      const [sx,sy,sz] = pocketInfo.size;
      v.addBox({ center:{x,y,z}, dimensions:{w:sx,h:sy,d:sz}, color:'red', alpha:0.2 });
      if (labelsVisible)
        v.addLabel('Binding Pocket', { position:{x,y,z}, fontColor:'black',
          backgroundColor:'white', backgroundOpacity:0.7 });
    }
    v.render();
  }, [receptorStyle, receptorColor, ligandStyle, ligandColor, surfaceVisible,
      pocketsVisible, labelsVisible, pocketInfo, moleculeType, isResultsMode]);

  // ── Load 3Dmol.js ────────────────────────────────────────────────────────
  const load3Dmol = useCallback(() => new Promise<void>((res, rej) => {
    if (window.$3Dmol) { res(); return; }
    const s = document.createElement('script');
    s.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
    s.onload = () => res(); s.onerror = () => rej(new Error('Failed to load 3Dmol.js'));
    document.head.appendChild(s);
  }), []);

  // ── Main render function ─────────────────────────────────────────────────
  const renderMolecule = useCallback(async (pathToLoad: string | null, extraPath?: string | null) => {
    if (!pathToLoad || !viewerRef.current) return;
    setIsLoading(true); setError(null); modelsLoaded.current = false;

    try {
      await load3Dmol();
      await new Promise(r => setTimeout(r, 0)); // flush DOM

      if (!viewerRef.current) throw new Error('Viewer detached');
      if (viewerInst.current) { viewerInst.current.clear(); }
      viewerInst.current = window.$3Dmol!.createViewer(viewerRef.current, { backgroundColor: 'white' });

      if (extraPath) {
        // Results mode: protein + docked ligand
        const [pData, lData] = await Promise.all([
          resolveStructureData(pathToLoad),
          resolveStructureData(extraPath),
        ]);
        viewerInst.current.addModel(pData.data, pData.format);
        viewerInst.current.addModel(lData.data, lData.format);
        setDataSource('Docking Result');
      } else {
        const { data, format, source } = await resolveStructureData(pathToLoad);
        viewerInst.current.addModel(data, format);
        setDataSource(source);
      }

      modelsLoaded.current = true;
      applyStyles();
      viewerInst.current.zoomTo();
      setIsLoading(false);
    } catch (e: any) {
      setError(e?.message || 'Visualization error');
      setIsLoading(false);
    }
  }, [load3Dmol, resolveStructureData, applyStyles]);

  // ── Effect: load when paths change ───────────────────────────────────────
  useEffect(() => {
    if (isResultsMode) {
      if (proteinPath && dockedPath) renderMolecule(proteinPath, dockedPath);
    } else {
      if (moleculePath) renderMolecule(moleculePath);
    }
    return () => {
      if (viewerInst.current) { viewerInst.current.clear(); viewerInst.current = null; }
    };
  }, [proteinPath, dockedPath, moleculePath, isResultsMode]);

  // ── Effect: re-apply styles when style toggles change ────────────────────
  useEffect(() => { applyStyles(); }, [applyStyles]);

  // ── Sync ligandList prop → currentLig ────────────────────────────────────
  useEffect(() => {
    if (ligandList?.length) setCurrentLig(c => c || ligandList[0].path);
    else setCurrentLig(undefined);
  }, [ligandList]);

  // ── Switch displayed ligand when selector changes ─────────────────────────
  const handleLigandSwitch = useCallback(async (path: string) => {
    setCurrentLig(path);
    onLigandChange?.(path);
    await renderMolecule(path);
  }, [renderMolecule, onLigandChange]);

  // ── Helpers ───────────────────────────────────────────────────────────────
  const infoChip = isResultsMode
    ? `Pose ${selectedPose + 1}${poseData?.[selectedPose]?.affinity ? ` · ${poseData[selectedPose].affinity} kcal/mol` : ''}`
    : dataSource;

  return (
    <div ref={containerRef} className="relative w-full rounded-lg border border-gray-200 bg-white overflow-hidden" style={{ height }}>

      {/* ── Ligand list selector (top-left, ligand-only views) ─────────────── */}
      {ligandList && ligandList.length > 1 && moleculeType === 'ligand' && (
        <div className="absolute top-2 left-2 z-10 flex items-center gap-2 bg-white/90 backdrop-blur border border-gray-200 rounded-lg px-2 py-1 shadow-sm">
          <span className="text-xs text-gray-500 font-medium">Ligand</span>
          <select
            className="text-xs border border-gray-300 rounded px-1 py-0.5 bg-white"
            value={currentLig}
            onChange={e => handleLigandSwitch(e.target.value)}
          >
            {ligandList.map(l => (
              <option key={l.path} value={l.path}>
                {l.name || l.path.split('/').pop() || l.path}
              </option>
            ))}
          </select>
          {onProceed && (
            <div className="flex gap-1 ml-1">
              <button onClick={() => onProceed('current', currentLig)} disabled={!currentLig}
                className="px-2 py-0.5 text-xs bg-blue-600 text-white rounded hover:bg-blue-700 disabled:opacity-40">
                This one
              </button>
              <button onClick={() => onProceed('all')}
                className="px-2 py-0.5 text-xs bg-green-600 text-white rounded hover:bg-green-700">
                All
              </button>
            </div>
          )}
        </div>
      )}

      {/* ── 3Dmol viewport ───────────────────────────────────────────────── */}
      <div ref={viewerRef} className="w-full h-full" />

      {/* ── Loading / error overlay ───────────────────────────────────────── */}
      {(isLoading || error) && (
        <div className="absolute inset-0 bg-white/80 backdrop-blur-sm flex flex-col items-center justify-center gap-3 p-4 text-center">
          {isLoading && (
            <>
              <svg className="w-8 h-8 animate-spin text-blue-500" fill="none" viewBox="0 0 24 24">
                <circle className="opacity-25" cx="12" cy="12" r="10" stroke="currentColor" strokeWidth="4"/>
                <path className="opacity-75" fill="currentColor" d="M4 12a8 8 0 018-8v8z"/>
              </svg>
              <p className="text-sm text-gray-600">Loading structure…</p>
            </>
          )}
          {error && (
            <>
              <span className="text-3xl">⚠️</span>
              <p className="text-sm font-medium text-red-600 max-w-xs">{error}</p>
              {moleculeType === 'ligand' && !moleculePath?.includes('/') && (
                <p className="text-xs text-gray-500">
                  Tip: click "Use Ligand ID" to confirm your compound before preview loads.
                </p>
              )}
            </>
          )}
        </div>
      )}

      {/* ── Toolbar (top-right) ───────────────────────────────────────────── */}
      {!isLoading && !error && (
        <>
          <div className="absolute top-2 right-2 flex items-center gap-1 p-1 bg-white/80 backdrop-blur-sm border border-gray-200 rounded-lg shadow-sm">
            {/* Style menu */}
            <div className="relative">
              <Btn title="Styles" onClick={() => setOpenMenu(o => o === 'style' ? null : 'style')}><StyleSVG/></Btn>
              {openMenu === 'style' && (
                <div className="absolute top-full right-0 mt-1 w-56 bg-white border rounded-lg shadow-lg p-3 z-20 text-sm space-y-2">
                  {[['Receptor', receptorStyle, setReceptorStyle, [['cartoon','Cartoon'],['sticks','Sticks'],['lines','Lines']]],
                    ['Receptor color', receptorColor, setReceptorColor, [['lightblue','Blue'],['spectrum','Spectrum'],['white','White'],['grey','Grey']]],
                    ['Ligand', ligandStyle, setLigandStyle, [['ballstick','Ball & Stick'],['stick','Stick'],['line','Line']]],
                    ['Ligand color', ligandColor, setLigandColor, [['Jmol','Jmol'],['cpk','CPK'],['grey','Grey']]],
                  ].map(([label, val, setter, opts]: any) => (
                    <div key={label} className="flex items-center justify-between gap-2">
                      <span className="text-gray-600 shrink-0">{label}</span>
                      <select value={val} onChange={e => setter(e.target.value)}
                        className="text-xs border border-gray-300 rounded px-1 py-0.5 bg-white flex-1 min-w-0">
                        {opts.map(([v,l]: string[]) => <option key={v} value={v}>{l}</option>)}
                      </select>
                    </div>
                  ))}
                </div>
              )}
            </div>

            {/* Options menu */}
            <div className="relative">
              <Btn title="Options" onClick={() => setOpenMenu(o => o === 'opts' ? null : 'opts')}><OptSVG/></Btn>
              {openMenu === 'opts' && (
                <div className="absolute top-full right-0 mt-1 w-44 bg-white border rounded-lg shadow-lg p-1 z-20 text-sm">
                  {[
                    ['Center view',         () => viewerInst.current?.zoomTo()],
                    [isSpinning ? 'Stop spin' : 'Spin', () => { setIsSpinning(s => !s); viewerInst.current?.spin(!isSpinning); }],
                    null,
                    [`${surfaceVisible?'Hide':'Show'} surface`, () => setSurfaceVisible(v => !v)],
                    [`${pocketsVisible?'Hide':'Show'} pocket`,  () => setPocketsVisible(v => !v)],
                    [`${labelsVisible ?'Hide':'Show'} labels`,  () => setLabelsVisible(v => !v)],
                    null,
                    ['Fullscreen', () => {
                      const el = containerRef.current;
                      if (!el) return;
                      document.fullscreenElement ? document.exitFullscreen().catch(()=>{}) : el.requestFullscreen().catch(()=>{});
                    }],
                    ['Download PNG', () => {
                      try {
                        const c = viewerRef.current?.querySelector('canvas') as HTMLCanvasElement;
                        if (!c) return;
                        const a = document.createElement('a'); a.href = c.toDataURL('image/png'); a.download = 'molecule.png'; a.click();
                      } catch {}
                    }],
                  ].map((item, i) => item === null
                    ? <hr key={i} className="my-1"/>
                    : <button key={i} onClick={item[1] as any} className="w-full text-left px-2 py-1 rounded hover:bg-gray-100">{item[0] as string}</button>
                  )}
                </div>
              )}
            </div>
          </div>

          {/* ── Info chip (bottom-left) ────────────────────────────────────── */}
          {infoChip && (
            <div className="absolute bottom-2 left-2 px-2 py-1 text-xs bg-white/80 backdrop-blur border border-gray-200 rounded-md shadow-sm">
              {infoChip}
            </div>
          )}
        </>
      )}
    </div>
  );
};