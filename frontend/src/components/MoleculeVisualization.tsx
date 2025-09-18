// src/components/MoleculeVisualization.tsx
import React, { useEffect, useRef, useState, useCallback } from 'react';

// --- Helper Components & Icons (for the UI) ---
const IconWrapper: React.FC<{ children: React.ReactNode; onClick?: () => void; title?: string }> = ({ children, onClick, title }) => (
    <button
        onClick={onClick}
        title={title}
        className="h-8 w-8 flex items-center justify-center rounded-md text-gray-600 hover:bg-gray-200 hover:text-gray-800 transition-colors"
    >
        {children}
    </button>
);
const StyleIcon = () => <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor"><path d="M17.414 2.586a2 2 0 00-2.828 0L7 10.172V13h2.828l7.586-7.586a2 2 0 000-2.828z" /><path fillRule="evenodd" d="M2 6a2 2 0 012-2h4a1 1 0 010 2H4v10h10v-4a1 1 0 112 0v4a2 2 0 01-2 2H4a2 2 0 01-2-2V6z" clipRule="evenodd" /></svg>;
const OptionsIcon = () => <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5" viewBox="0 0 20 20" fill="currentColor"><path fillRule="evenodd" d="M11.49 3.17c-.38-1.56-2.6-1.56-2.98 0a1.532 1.532 0 01-2.286.948c-1.372-.836-2.942.734-2.106 2.106.54.886.061 2.042-.947 2.287-1.561.379-1.561 2.6 0 2.978a1.532 1.532 0 01.947 2.287c-.836 1.372.734 2.942 2.106 2.106a1.532 1.532 0 012.287.947c.379 1.561 2.6 1.561 2.978 0a1.532 1.532 0 012.287-.947c1.372.836 2.942-.734 2.106-2.106a1.532 1.532 0 01-.947-2.287c1.561-.379 1.561-2.6 0-2.978a1.532 1.532 0 01-.947-2.287c.836-1.372-.734-2.942-2.106-2.106a1.532 1.532 0 01-2.287-.947zM10 13a3 3 0 100-6 3 3 0 000 6z" clipRule="evenodd" /></svg>;

// --- Type Definitions ---
declare global { interface Window { $3Dmol?: any } }
type MolFormat = 'pdb' | 'pdbqt' | 'sdf' | 'mol2';
interface Pocket { center: number[]; size: number[]; }
interface PoseRow { affinity: string; }

interface MoleculeVisualizationProps {
    moleculePath?: string | null;
    moleculeType?: 'protein' | 'ligand';
    proteinPath?: string | null;
    dockedPath?: string | null;
    height?: string;
    pocketInfo?: Pocket;
    showPockets?: boolean;
    selectedPose?: number;
    isResultsMode?: boolean;
    poseData?: PoseRow[];
}

// --- Main Component ---
export const MoleculeVisualization: React.FC<MoleculeVisualizationProps> = ({
    moleculePath,
    moleculeType = 'protein',
    proteinPath,
    dockedPath,
    height = '500px',
    pocketInfo,
    showPockets = true,
    isResultsMode = false,
    selectedPose = 0,
    poseData,
}) => {
    const containerRef = useRef<HTMLDivElement | null>(null);
    const viewerRef = useRef<HTMLDivElement | null>(null);
    const viewer = useRef<any>(null);
    const modelsLoaded = useRef<boolean>(false); // <-- Key change: Track model state

    // State
    const [isLoading, setIsLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [dataSource, setDataSource] = useState('');
    const [openMenu, setOpenMenu] = useState<string | null>(null);

    // Style & Feature State
    const [surfaceVisible, setSurfaceVisible] = useState(true);
    const [labelsVisible, setLabelsVisible] = useState(false);
    const [pocketsVisible, setPocketsVisible] = useState(showPockets);
    const [isSpinning, setIsSpinning] = useState(false);
    const [ligandStyle, setLigandStyle] = useState<'ballstick' | 'stick' | 'line'>('ballstick');
    const [receptorStyle, setReceptorStyle] = useState<'cartoon' | 'sticks' | 'lines'>('cartoon');
    const [receptorColor, setReceptorColor] = useState<'lightblue' | 'spectrum' | 'white' | 'grey'>('lightblue');
    const [ligandColor, setLigandColor] = useState<'Jmol' | 'cpk' | 'grey'>('Jmol');

    const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

    const resolveStructureData = useCallback(async (input: string): Promise<{ data: string; format: string; source: string }> => {
        // This function remains the same as it's robust
        if (!input) throw new Error('No input provided');
        const trimmedInput = input.trim();
        if (trimmedInput.includes('/') || trimmedInput.includes('\\')) {
            const url = trimmedInput.startsWith('/api/download') ? `${apiBase}${trimmedInput}` : `${apiBase}/api/download?path=${encodeURIComponent(trimmedInput)}`;
            const res = await fetch(url);
            if (!res.ok) throw new Error(`Download failed: ${res.status}`);
            const text = await res.text();
            const ext = (trimmedInput.split('.').pop() || '').toLowerCase();
            return { data: text, format: ['pdb', 'pdbqt', 'sdf', 'mol2'].includes(ext) ? ext as MolFormat : 'pdb', source: 'Server File' };
        }
        if (/^[1-9][A-Za-z0-9]{3}$/i.test(trimmedInput)) {
            const pdbUrl = `https://files.rcsb.org/download/${trimmedInput.toUpperCase()}.pdb`;
            const res = await fetch(pdbUrl);
            if (res.ok) return { data: await res.text(), format: 'pdb', source: 'RCSB PDB' };
            throw new Error('Failed to fetch from RCSB PDB');
        }
        const isNumeric = /^[0-9]+$/.test(trimmedInput);
        const searchType = isNumeric ? 'cid' : 'name';
        const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${searchType}/${encodeURIComponent(trimmedInput)}/SDF`;
        const proxyUrl = `https://api.allorigins.win/raw?url=${encodeURIComponent(pubchemUrl)}`;
        const res = await fetch(proxyUrl);
        if (res.ok) {
             const data = await res.text();
             if (data && !data.includes('Status: 404')) return { data, format: 'sdf', source: 'PubChem' };
        }
        throw new Error(`Could not resolve: "${trimmedInput}"`);
    }, [apiBase]);

    const applyStylesAndExtras = useCallback(() => {
        if (!viewer.current || !modelsLoaded.current) return; // <-- Key change: Guard clause

        // Receptor Style
        const rColorSpec = { color: receptorColor === 'spectrum' ? 'spectrum' : receptorColor };
        if (receptorStyle === 'cartoon') viewer.current.setStyle({ hetflag: false }, { cartoon: { ...rColorSpec, opacity: 0.9 } });
        else if (receptorStyle === 'sticks') viewer.current.setStyle({ hetflag: false }, { stick: { radius: 0.2, ...rColorSpec } });
        else viewer.current.setStyle({ hetflag: false }, { line: rColorSpec });

        // Ligand Style
        const scheme = ligandColor === 'grey' ? 'greyCarbon' : ligandColor;
        const target = isResultsMode ? { model: 1 } : {};
        if (ligandStyle === 'ballstick') viewer.current.setStyle(target, { stick: { radius: 0.3, colorscheme: scheme }, sphere: { radius: 0.5, colorscheme: scheme } });
        else if (ligandStyle === 'stick') viewer.current.setStyle(target, { stick: { radius: 0.25, colorscheme: scheme } });
        else viewer.current.setStyle(target, { line: { colorscheme: scheme } });

        // Surfaces
        viewer.current.removeAllSurfaces();
        if (surfaceVisible && (moleculeType === 'protein' || isResultsMode)) {
            viewer.current.addSurface('MS', { opacity: 0.6, color: 'white' }, { model: 0 });
        }
        
        // Pockets & Labels
        viewer.current.removeAllShapes();
        viewer.current.removeAllLabels();
        if (pocketsVisible && pocketInfo) {
            const [x, y, z] = pocketInfo.center;
            const [sx, sy, sz] = pocketInfo.size;
            viewer.current.addBox({ center: { x, y, z }, dimensions: { w: sx, h: sy, d: sz }, color: 'red', alpha: 0.2 });
            if (labelsVisible) viewer.current.addLabel('Binding Pocket', { position: { x, y, z }, fontColor: 'black', backgroundColor: 'white', backgroundOpacity: 0.7 });
        }
        
        viewer.current.render();
    }, [isResultsMode, receptorStyle, receptorColor, ligandStyle, ligandColor, surfaceVisible, pocketsVisible, labelsVisible, pocketInfo, moleculeType]);

    // Initial load and reload on path changes
    useEffect(() => {
        // Encapsulate all 3Dmol logic in one effect
        const initAndRender = async () => {
            if (!containerRef.current || !(isResultsMode ? proteinPath && dockedPath : moleculePath)) return;

            setIsLoading(true);
            setError(null);
            modelsLoaded.current = false;

            try {
                // Ensure 3Dmol script is loaded
                if (!window.$3Dmol) {
                    await new Promise<void>((resolve, reject) => {
                        const s = document.createElement('script');
                        s.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
                        s.onload = () => resolve();
                        s.onerror = () => reject(new Error('Failed to load 3Dmol.js'));
                        document.head.appendChild(s);
                    });
                }
                
                // Key change: Use setTimeout to ensure the container div has its dimensions
                setTimeout(async () => {
                    try {
                        if (!viewerRef.current) throw new Error("Viewer element disappeared");
                        viewer.current = window.$3Dmol.createViewer(viewerRef.current, { backgroundColor: 'white' });

                        if (isResultsMode) {
                            if (!proteinPath || !dockedPath) throw new Error('Protein or Ligand path missing.');
                            const [pData, lData] = await Promise.all([ resolveStructureData(proteinPath), resolveStructureData(dockedPath) ]);
                            viewer.current.addModel(pData.data, pData.format);
                            viewer.current.addModel(lData.data, lData.format);
                            setDataSource('Docking Result');
                        } else {
                            if (!moleculePath) throw new Error('Molecule path missing.');
                            const { data, format, source } = await resolveStructureData(moleculePath);
                            viewer.current.addModel(data, format);
                            setDataSource(source);
                        }

                        modelsLoaded.current = true; // <-- Signal that models are ready
                        applyStylesAndExtras(); // <-- Now apply styles
                        viewer.current.zoomTo();
                        setIsLoading(false);
                    } catch (e: any) {
                        setError(e.message || 'Error during rendering.');
                        setIsLoading(false);
                    }
                }, 0); // <-- Defer execution

            } catch (e: any) {
                setError(e.message);
                setIsLoading(false);
            }
        };

        initAndRender();
        
        // Cleanup on unmount
        return () => {
            if (viewer.current) {
                viewer.current.clear();
                viewer.current = null;
            }
        };
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [proteinPath, dockedPath, moleculePath, isResultsMode]); // Re-run only when paths change

    // Effect for re-applying styles when UI controls change
    useEffect(() => {
        applyStylesAndExtras();
    }, [applyStylesAndExtras]);

    // UI Handlers... (no changes here)
    const toggleFullscreen = () => { /* ... */ };
    const handleDownloadImage = () => { /* ... */ };

    return (
        <div ref={containerRef} className="relative w-full rounded-lg border border-gray-200 bg-white overflow-hidden" style={{ height }}>
           {/* JSX for UI remains the same */}
            <div ref={viewerRef} className="w-full h-full" />
            
            {(isLoading || error) && (
                <div className="absolute inset-0 bg-white/70 backdrop-blur-sm flex items-center justify-center text-center p-4">
                    {isLoading && <p>Loading Visualization...</p>}
                    {error && <p className="font-medium text-red-600">Error: {error}</p>}
                </div>
            )}

            {!isLoading && !error && (
                <>
                    {/* Toolbar */}
                    <div className="absolute top-2 right-2 flex items-center gap-1 p-1 bg-white/80 backdrop-blur-sm border border-gray-200 rounded-lg shadow-sm">
                        {/* Style Menu */}
                        <div>
                            <IconWrapper onClick={() => setOpenMenu(o => o === 'style' ? null : 'style')} title="Styles"><StyleIcon /></IconWrapper>
                            {openMenu === 'style' && (
                                <div className="absolute top-full right-0 mt-2 w-60 bg-white border rounded-md shadow-lg p-3 z-10 text-sm space-y-3">
                                    <div className="flex items-center justify-between">
                                        <label>Receptor</label>
                                        <select value={receptorStyle} onChange={e => setReceptorStyle(e.target.value as any)} className="text-xs border-gray-300 rounded-md ml-2">
                                            <option value="cartoon">Cartoon</option><option value="sticks">Sticks</option><option value="lines">Lines</option>
                                        </select>
                                    </div>
                                    <div className="flex items-center justify-between">
                                        <label>Receptor Color</label>
                                        <select value={receptorColor} onChange={e => setReceptorColor(e.target.value as any)} className="text-xs border-gray-300 rounded-md ml-2">
                                            <option value="lightblue">Blue</option><option value="spectrum">Spectrum</option><option value="white">White</option><option value="grey">Grey</option>
                                        </select>
                                    </div>
                                    <hr/>
                                    <div className="flex items-center justify-between">
                                        <label>Ligand</label>
                                        <select value={ligandStyle} onChange={e => setLigandStyle(e.target.value as any)} className="text-xs border-gray-300 rounded-md ml-2">
                                            <option value="ballstick">Ball & Stick</option><option value="stick">Stick</option><option value="line">Line</option>
                                        </select>
                                    </div>
                                    <div className="flex items-center justify-between">
                                        <label>Ligand Color</label>
                                        <select value={ligandColor} onChange={e => setLigandColor(e.target.value as any)} className="text-xs border-gray-300 rounded-md ml-2">
                                            <option value="Jmol">Jmol</option><option value="cpk">CPK</option><option value="grey">Grey</option>
                                        </select>
                                    </div>
                                </div>
                            )}
                        </div>
                        {/* Options Menu */}
                        <div>
                            <IconWrapper onClick={() => setOpenMenu(o => o === 'options' ? null : 'options')} title="Options"><OptionsIcon /></IconWrapper>
                            {openMenu === 'options' && (
                                <div className="absolute top-full right-0 mt-2 w-48 bg-white border rounded-md shadow-lg p-2 z-10 text-sm">
                                    <button onClick={() => viewer.current?.zoomTo()} className="w-full text-left p-1 rounded hover:bg-gray-100">Center View</button>
                                    <button onClick={() => { setIsSpinning(s => !s); viewer.current?.spin(!isSpinning); }} className="w-full text-left p-1 rounded hover:bg-gray-100">{isSpinning ? 'Stop Spin' : 'Spin'}</button>
                                    <hr className="my-1"/>
                                    <button onClick={() => setSurfaceVisible(v => !v)} className="w-full text-left p-1 rounded hover:bg-gray-100">{surfaceVisible ? 'Hide' : 'Show'} Surface</button>
                                    <button onClick={() => setPocketsVisible(v => !v)} className="w-full text-left p-1 rounded hover:bg-gray-100">{pocketsVisible ? 'Hide' : 'Show'} Pocket</button>
                                    <button onClick={() => setLabelsVisible(v => !v)} className="w-full text-left p-1 rounded hover:bg-gray-100">{labelsVisible ? 'Hide' : 'Show'} Labels</button>
                                    <hr className="my-1"/>
                                    <button onClick={toggleFullscreen} className="w-full text-left p-1 rounded hover:bg-gray-100">Fullscreen</button>
                                    <button onClick={handleDownloadImage} className="w-full text-left p-1 rounded hover:bg-gray-100">Download PNG</button>
                                </div>
                            )}
                        </div>
                    </div>
                    {/* Info Chip */}
                    <div className="absolute bottom-2 left-2 px-2 py-1 text-xs bg-white/80 backdrop-blur rounded-md border border-gray-200 shadow-sm">
                        {isResultsMode ? `Pose ${selectedPose + 1}${poseData && poseData[selectedPose]?.affinity ? ` • ${poseData[selectedPose].affinity} kcal/mol` : ''}` : dataSource}
                    </div>
                </>
            )}
        </div>
    );
};