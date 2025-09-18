// src/components/MoleculeVisualization.tsx
import React, { useEffect, useRef, useState, useCallback } from 'react';

// Safe 3Dmol types
declare global {
  interface Window {
    $3Dmol?: {
      createViewer: (element: HTMLElement, config?: Record<string, unknown>) => {
        setBackgroundColor: (color: string) => void;
        clear: () => void;
        addModel: (data: string, format: string) => any;
        setStyle: (selector: Record<string, unknown>, style: Record<string, unknown>) => void;
        addSphere: (sphere: { center: { x: number; y: number; z: number }; radius: number; color: string; alpha?: number }) => void;
        addBox: (box: { center: { x: number; y: number; z: number }; dimensions: { w: number; h: number; d: number }; color?: string; alpha?: number; wireframe?: boolean }) => void;
        zoomTo: (selector?: Record<string, unknown>, animationDuration?: number, fixedPath?: boolean) => void;
        render: () => void;
        resize: () => void;
        spin: (enable: boolean) => void;
        zoom: (factor: number, animationDuration?: number) => void;
        slab?: (val: number) => void;
        setSlab?: (near: number, far: number) => void;
      } | null;
      rasmolElementColors: Record<string, string>;
    };
  }
}

interface PoseRow {
  pose: number;
  affinity: string;
  rmsd_lb: string;
  rmsd_ub: string;
}

interface Pocket {
  center: number[];
  size: number[];
  confidence: string;
  method: string;
}

interface MoleculeVisualizationProps {
  // Generic path (non-results mode)
  moleculePath?: string | null;

  // Results mode: pass both files
  proteinPath?: string | null;     // prepared protein PDB/PDBQT
  dockedPath?: string | null;      // docked poses PDBQT

  moleculeType?: 'protein' | 'ligand' | 'complex';
  color?: string;
  height?: string;

  pocketInfo?: Pocket;
  allPockets?: Pocket[];

  selectedPose?: number;
  showPockets?: boolean;
  jobId?: string;
  isResultsMode?: boolean;
  poseData?: PoseRow[];
}

export const MoleculeVisualization: React.FC<MoleculeVisualizationProps> = ({
  moleculePath,
  proteinPath,
  dockedPath,
  moleculeType = 'protein',
  color = '#6b7280',
  height = '400px',
  pocketInfo,
  allPockets,
  selectedPose = 0,
  showPockets = true,
  jobId,
  isResultsMode = false,
  poseData
}) => {
  const viewerRef = useRef<HTMLDivElement | null>(null);
  const containerRef = useRef<HTMLDivElement | null>(null);
  const viewerInstanceRef = useRef<any>(null);
  const currentFormatRef = useRef<string>('pdb');

  const [isLoading, setIsLoading] = useState(false);
  const [hasViewer, setHasViewer] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [dataSource, setDataSource] = useState('');
  const [currentData, setCurrentData] = useState('');
  const [pocketsVisible, setPocketsVisible] = useState(showPockets);

  const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

  const isServerFilePath = (p: string) => {
    if (!p) return false;
    const s = p.trim();
    return s.startsWith('/api/download')
      || s.includes('temp_runs')
      || s.includes('data/results')
      || s.includes('docking')
      || s.includes('uploads')
      || s.includes('prepared')
      || s.includes('data/proteins')
      || s.includes('data/ligands');
  }; // [file:1]

  const buildDownloadUrl = (p: string) => {
    const s = p.trim();
    if (s.startsWith('/api/download')) return `${apiBase}${s}`;
    return `${apiBase}/api/download?path=${encodeURIComponent(s)}`;
  }; // [file:1]

  const sniffFormat = (filePath: string, text: string): string => {
    const ext = (filePath.split('.').pop() || '').toLowerCase();
    if (ext === 'pdb') return 'pdb';
    if (ext === 'pdbqt') return 'pdbqt';
    if (ext === 'sdf') return 'sdf';
    if (ext === 'mol2') return 'mol2';
    const t = text || '';
    if (t.includes('@<TRIPOS>MOLECULE')) return 'mol2';
    if (t.includes('M  END')) return 'sdf';
    if (t.includes('REMARK') && t.includes('VINA')) return 'pdbqt';
    if (t.includes('ATOM') || t.includes('HETATM')) return 'pdb';
    return 'pdb';
  }; // [file:1]

  const fetchServerFile = async (filePath: string): Promise<{ data: string; format: string; source: string }> => {
    const url = buildDownloadUrl(filePath);
    const res = await fetch(url);
    if (!res.ok) throw new Error(`Download failed ${res.status}`);
    const text = await res.text();
    const format = sniffFormat(filePath, text);
    return { data: text, format, source: 'Server File' };
  }; // [file:1]

  const resolveStructureData = async (input: string): Promise<{ data: string; format: string; source: string }> => {
    if (!input || input.trim() === '') throw new Error('Empty input provided');

    const trimmedInput = input.trim();

    // Always route server-side files through API
    if (isServerFilePath(trimmedInput)) {
      const r = await fetchServerFile(trimmedInput);
      return { ...r, source: isResultsMode ? 'Docking Results' : r.source };
    }

    // If looks like URL, fetch directly
    if (/^https?:\/\//i.test(trimmedInput)) {
      const res = await fetch(trimmedInput);
      if (!res.ok) throw new Error(`URL fetch failed: ${res.status}`);
      const text = await res.text();
      const format = sniffFormat(trimmedInput, text);
      return { data: text, format, source: 'URL' };
    }

    // Raw PDB text
    if (trimmedInput.includes('ATOM') || trimmedInput.includes('HETATM')) {
      return { data: trimmedInput, format: 'pdb', source: 'Direct PDB Data' };
    }

    // PDB ID
    if (/^[1-9][A-Za-z0-9]{3}$/i.test(trimmedInput)) {
      const pdbUrl = `https://files.rcsb.org/download/${trimmedInput.toUpperCase()}.pdb`;
      const response = await fetch(`https://api.allorigins.win/raw?url=${encodeURIComponent(pdbUrl)}`);
      if (response.ok) {
        const data = await response.text();
        if (data.includes('ATOM') || data.includes('HETATM')) {
          return { data, format: 'pdb', source: 'RCSB PDB Database' };
        }
      }
      throw new Error('Failed to resolve PDB ID');
    }

    // PubChem name or CID -> SDF
    {
      const isNumeric = /^[0-9]+$/.test(trimmedInput);
      const searchType = isNumeric ? 'cid' : 'name';
      const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${searchType}/${encodeURIComponent(trimmedInput)}/SDF`;
      const response = await fetch(`https://api.allorigins.win/raw?url=${encodeURIComponent(pubchemUrl)}`);
      if (response.ok) {
        const data = await response.text();
        if (data.includes('M  END')) {
          return { data, format: 'sdf', source: 'PubChem Database' };
        }
      }
    }

    throw new Error(`Could not resolve structure: "${trimmedInput}"`);
  }; // [file:1]

  const extractPoseFromData = (data: string, poseIndex: number): string => {
    if (!isResultsMode || poseIndex === 0) return data;
    const fmt = currentFormatRef.current;
    if ((fmt === 'pdb' || fmt === 'pdbqt') && data.includes('MODEL')) {
      const chunks = data.split(/MODEL\s+\d+/);
      if (chunks.length > poseIndex + 1) {
        const selected = chunks[poseIndex + 1];
        const endIndex = selected.indexOf('ENDMDL');
        if (endIndex !== -1) {
          return `MODEL ${poseIndex + 1}\n${selected.substring(0, endIndex)}\nENDMDL\nEND`;
        }
      }
    }
    return data;
  }; // [file:1]

  const addPocketVisualization = useCallback((viewer: any) => {
    if (!pocketsVisible || !viewer) return;
    if (pocketInfo) {
      const [x, y, z] = pocketInfo.center;
      const [sx, sy, sz] = pocketInfo.size;
      viewer.addBox({
        center: { x, y, z },
        dimensions: { w: sx, h: sy, d: sz },
        color: '#ef4444',
        alpha: 0.25,
        wireframe: true
      });
      viewer.addSphere({ center: { x, y, z }, radius: 2.5, color: '#ef4444', alpha: 0.7 });
    }
    if (allPockets && allPockets.length > 1) {
      const colors = ['#2563eb', '#10b981', '#f59e0b', '#8b5cf6', '#eab308', '#06b6d4'];
      allPockets.slice(1).forEach((p, idx) => {
        const [x, y, z] = p.center;
        const [sx, sy, sz] = p.size;
        const c = colors[idx % colors.length];
        viewer.addBox({
          center: { x, y, z },
          dimensions: { w: sx, h: sy, d: sz },
          color: c,
          alpha: 0.2,
          wireframe: true
        });
        viewer.addSphere({ center: { x, y, z }, radius: 2.0, color: c, alpha: 0.5 });
      });
    }
  }, [pocketInfo, allPockets, pocketsVisible]); // [file:1]

  const load3DmolIfNeeded = async () => {
    if (window.$3Dmol) return;
    await new Promise<void>((resolve, reject) => {
      const script = document.createElement('script');
      script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
      script.onload = () => setTimeout(() => window.$3Dmol ? resolve() : reject(new Error('3Dmol not available')), 150);
      script.onerror = () => reject(new Error('Failed to load 3Dmol script'));
      document.head.appendChild(script);
    });
  }; // [file:1]

  const loadTextViaApi = async (path: string) => {
    const url = path.startsWith('/api/download') ? `${apiBase}${path}` : `${apiBase}/api/download?path=${encodeURIComponent(path)}`;
    const r = await fetch(url);
    if (!r.ok) throw new Error(`Fetch failed: ${r.status}`);
    return await r.text();
  }; // [file:1]

  const loadResultsComplex = async (viewer: any) => {
    if (!proteinPath || !dockedPath) throw new Error('Missing protein or docked file');
    // Fetch both texts
    const [protText, dockText] = await Promise.all([loadTextViaApi(proteinPath), loadTextViaApi(dockedPath)]);
    const protFmt = sniffFormat(proteinPath, protText);
    const dockFmt = 'pdbqt';
    currentFormatRef.current = dockFmt;
    setCurrentData(dockText);

    const poseText = extractPoseFromData(dockText, selectedPose);

    // Build scene
    viewer.clear();
    const protModel = viewer.addModel(protText, protFmt);
    if (!protModel) throw new Error('Protein model failed to load');

    // Protein style: cartoon
    viewer.setStyle({}, { cartoon: { color: 'lightblue', opacity: 0.95 } });

    const ligModel = viewer.addModel(poseText, dockFmt);
    if (!ligModel) throw new Error('Ligand model failed to load');

    // Ligand style: ball-and-stick
    viewer.setStyle({ model: ligModel }, {
      stick: { radius: 0.65, colorscheme: 'Jmol' },
      sphere: { radius: 0.7, colorscheme: 'Jmol', opacity: 0.9 }
    });

    // Highlight pocket region and nearby residues
    if (pocketInfo) {
      const [x, y, z] = pocketInfo.center;
      const [sx, sy, sz] = pocketInfo.size;
      viewer.setStyle({
        and: [
          { atom: 'CA' },
          { x: { $gte: x - sx/2, $lte: x + sx/2 } },
          { y: { $gte: y - sy/2, $lte: y + sy/2 } },
          { z: { $gte: z - sz/2, $lte: z + sz/2 } }
        ]
      }, { cartoon: { color: 'orange', opacity: 1.0 } });
    }

    if (pocketsVisible) addPocketVisualization(viewer);

    // Camera framing similar to CB-Dock2
    viewer.setBackgroundColor('#ffffff');
    if (pocketInfo) {
      const [x, y, z] = pocketInfo.center;
      const [sx, sy, sz] = pocketInfo.size;
      viewer.zoomTo({
        and: [
          { x: { $gte: x - sx, $lte: x + sx } },
          { y: { $gte: y - sy, $lte: y + sy } },
          { z: { $gte: z - sz, $lte: z + sz } }
        ]
      }, 650);
    } else {
      viewer.zoomTo({}, 650);
    }
    // Clipping
    viewer.setSlab && viewer.setSlab(5, 140);
    viewer.zoom(0.92, 500);
    viewer.render();

    setHasViewer(true);
    setDataSource('Docking Results (protein+ligand)');
  }; // [file:1]

  // Generic, non-results loader
  const loadSingleStructure = async () => {
    if (!moleculePath) return;
    await load3DmolIfNeeded();
    if (!viewerRef.current || !containerRef.current) throw new Error('Viewer elements not ready');

    // Size
    const rect = containerRef.current.getBoundingClientRect();
    const w = Math.max(rect.width, 320);
    const h = Math.max(rect.height, 320);
    viewerRef.current.style.width = `${w}px`;
    viewerRef.current.style.height = `${h}px`;

    const viewer = window.$3Dmol!.createViewer(viewerRef.current!, { backgroundColor: 'white' });
    viewerInstanceRef.current = viewer;
    const { data, format, source } = await resolveStructureData(moleculePath);
    currentFormatRef.current = format;
    setDataSource(source);
    setCurrentData(data);

    viewer.clear();
    const model = viewer.addModel(data, format);
    if (!model) throw new Error('Failed to create model from data');

    if (format === 'pdb' && moleculeType === 'protein') {
      viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    } else {
      viewer.setStyle({}, { stick: { radius: 0.4, colorscheme: 'Jmol' }, sphere: { radius: 0.45, colorscheme: 'Jmol' } });
    }

    addPocketVisualization(viewer);

    viewer.zoomTo();
    viewer.zoom(1.1, 500);
    viewer.render();
    setHasViewer(true);
  }; // [file:1]

  // Initial load
  useEffect(() => {
    const init = async () => {
      if (!(moleculePath || (isResultsMode && proteinPath && dockedPath))) return;
      try {
        setIsLoading(true);
        setError(null);
        await load3DmolIfNeeded();

        // Prepare container size early
        if (viewerRef.current && containerRef.current) {
          const rect = containerRef.current.getBoundingClientRect();
          viewerRef.current.style.width = `${Math.max(rect.width, 320)}px`;
          viewerRef.current.style.height = `${Math.max(rect.height, 320)}px`;
        }

        const viewer = window.$3Dmol!.createViewer(viewerRef.current!, { backgroundColor: 'white' });
        viewerInstanceRef.current = viewer;

        if (isResultsMode && proteinPath && dockedPath) {
          await loadResultsComplex(viewer);
        } else if (moleculePath) {
          await loadSingleStructure();
        }
      } catch (e: any) {
        console.error('Viewer initialization failed:', e);
        setError(e?.message || 'Failed to load viewer');
        setHasViewer(false);
        setDataSource('');
      } finally {
        setIsLoading(false);
      }
    };
    init();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [moleculePath, proteinPath, dockedPath, moleculeType, isResultsMode]); // [file:1]

  // Pose change in results mode
  useEffect(() => {
    if (!viewerInstanceRef.current || !currentData || !isResultsMode || !dockedPath || !proteinPath) return;
    const rerender = async () => {
      try {
        const viewer = viewerInstanceRef.current;
        await loadResultsComplex(viewer);
      } catch (e) {
        console.warn('Pose switch failed', e);
      }
    };
    rerender();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selectedPose, pocketsVisible]); // [file:1]

  // Pocket toggle re-render (for non-results mode)
  useEffect(() => {
    if (!viewerInstanceRef.current || !currentData || isResultsMode) return;
    try {
      const viewer = viewerInstanceRef.current;
      viewer.clear();
      const model = viewer.addModel(currentData, currentFormatRef.current);
      if (model) {
        if (currentFormatRef.current === 'pdb' && moleculeType === 'protein') {
          viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
        } else {
          viewer.setStyle({}, { stick: { radius: 0.4, colorscheme: 'Jmol' }, sphere: { radius: 0.45, colorscheme: 'Jmol' } });
        }
        if (pocketsVisible) addPocketVisualization(viewer);
        viewer.render();
      }
    } catch (e) {
      console.warn('Pocket toggle re-render failed', e);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pocketsVisible]); // [file:1]

  const togglePockets = () => setPocketsVisible(!pocketsVisible); // [file:1]

  const isFilePath = (moleculePath || proteinPath || dockedPath || '').includes('/') || (moleculePath || proteinPath || dockedPath || '').includes('\\'); // [file:1]
  const displayName = isFilePath ? (isResultsMode ? 'Docking Results' : 'Uploaded File') : (moleculePath || ''); // [file:1]

  return (
    <div className="w-full">
      <div ref={containerRef} className="w-full relative rounded border border-gray-200 bg-white" style={{ height }}>
        <div ref={viewerRef} className="w-full h-full" />

        {/* Compact top-left toolbar */}
        {hasViewer && !isLoading && !error && (
          <div className="absolute top-2 left-2 flex items-center gap-2">
            <span className="px-2 py-1 text-xs bg-white/80 backdrop-blur rounded border border-gray-200 shadow-sm">
              {isResultsMode ? `Pose ${selectedPose + 1}${poseData && poseData[selectedPose]?.affinity ? ` • ${poseData[selectedPose].affinity} kcal/mol` : ''}` : '3D'}
            </span>
            {(pocketInfo || (allPockets && allPockets.length > 0)) && (
              <button className="px-2 py-1 text-xs bg-white/80 rounded border hover:bg-white" onClick={togglePockets}>
                {pocketsVisible ? 'Hide Pockets' : 'Show Pockets'}
              </button>
            )}
            <button className="px-2 py-1 text-xs bg-white/80 rounded border hover:bg-white" onClick={() => viewerInstanceRef.current?.spin?.(true)}>Rotate</button>
            <button className="px-2 py-1 text-xs bg-white/80 rounded border hover:bg-white" onClick={() => { viewerInstanceRef.current?.zoomTo?.(); viewerInstanceRef.current?.render?.(); }}>Reset</button>
          </div>
        )}

        {/* Minimal error/loading chip, bottom-left */}
        {(isLoading || error || (!hasViewer && !isLoading)) && (
          <div className="absolute bottom-2 left-2 px-2 py-1 text-xs bg-white/80 backdrop-blur rounded border border-gray-200 shadow-sm">
            {error ? `Error: ${error}` : isLoading ? `Loading ${isResultsMode ? 'results' : 'structure'}...` : (dataSource || displayName)}
          </div>
        )}
      </div>
    </div>
  );
};
