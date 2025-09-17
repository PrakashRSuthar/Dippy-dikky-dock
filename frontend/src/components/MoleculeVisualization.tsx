import React, { useEffect, useRef, useState, useCallback } from 'react';

// Safe 3Dmol types
declare global {
  interface Window {
    $3Dmol?: {
      createViewer: (element: HTMLElement, config?: Record<string, unknown>) => {
        setBackgroundColor: (color: string) => void;
        clear: () => void;
        addModel: (data: string, format: string) => unknown;
        setStyle: (selector: Record<string, unknown>, style: Record<string, unknown>) => void;
        addSphere: (sphere: { center: { x: number; y: number; z: number }; radius: number; color: string; alpha?: number }) => void;
        addBox: (box: { center: { x: number; y: number; z: number }; dimensions: { w: number; h: number; d: number }; color?: string; alpha?: number; wireframe?: boolean }) => void;
        zoomTo: (selector?: Record<string, unknown>, animationDuration?: number, fixedPath?: boolean) => void;
        render: () => void;
        resize: () => void;
        setStyle: (selector: Record<string, unknown>, style: Record<string, unknown>) => void;
        spin: (enable: boolean) => void;
        zoom: (factor: number, animationDuration?: number) => void;
      } | null;
      rasmolElementColors: Record<string, string>;
    };
  }
}

interface MoleculeVisualizationProps {
  moleculePath?: string | null;
  moleculeType?: 'protein' | 'ligand' | 'complex';
  color?: string;
  height?: string;
  // Docking-specific props
  pocketInfo?: {
    center: number[];
    size: number[];
    confidence: string;
    method: string;
  };
  allPockets?: Array<{
    center: number[];
    size: number[];
    confidence: string;
    method: string;
  }>;
  selectedPose?: number;
  showPockets?: boolean;
  jobId?: string;
  isResultsMode?: boolean;
  // NEW: Pose data for multi-conformer files
  poseData?: Array<{ pose: number; affinity: string; rmsd_lb: string; rmsd_ub: string }>;
}

export const MoleculeVisualization: React.FC<MoleculeVisualizationProps> = ({ 
  moleculePath, 
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
  const viewerRef = useRef<HTMLDivElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [hasViewer, setHasViewer] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [dataSource, setDataSource] = useState<string>('');
  const [currentData, setCurrentData] = useState<string>('');
  const [pocketsVisible, setPocketsVisible] = useState(showPockets);

  // Docking result file fetcher
  const fetchDockingResultFile = async (filePath: string): Promise<{data: string, format: string, source: string}> => {
    try {
      if (filePath.startsWith('/') || filePath.includes('temp_runs') || filePath.includes('data/results') || filePath.includes('docking')) {
        const apiBase = import.meta.env.VITE_API_BASE || 'http://localhost:8000';
        const downloadUrl = `${apiBase}/api/download?path=${encodeURIComponent(filePath)}`;
        
        const response = await fetch(downloadUrl);
        if (!response.ok) {
          throw new Error(`Failed to fetch docking result: ${response.status}`);
        }
        
        const data = await response.text();
        
        let format = 'pdb';
        if (filePath.endsWith('.pdbqt')) {
          format = 'pdbqt';
        } else if (filePath.endsWith('.sdf')) {
          format = 'sdf';
        }
        
        return { data, format, source: 'Docking Results' };
      }
      
      throw new Error('Invalid docking result file path');
    } catch (err) {
      console.error('Docking file fetch failed:', err);
      throw err;
    }
  };

  // Structure resolution
  const resolveStructureData = async (input: string): Promise<{data: string, format: string, source: string}> => {
    if (!input || input.trim() === '') {
      throw new Error('Empty input provided');
    }

    const trimmedInput = input.trim();
    
    // Docking result file path
    if (isResultsMode && (trimmedInput.includes('temp_runs') || trimmedInput.includes('data/results') || trimmedInput.includes('docking'))) {
      return await fetchDockingResultFile(trimmedInput);
    }

    // Regular file path
    if (trimmedInput.includes('/') || trimmedInput.includes('\\')) {
      try {
        const response = await fetch(trimmedInput);
        if (!response.ok) throw new Error(`File fetch failed: ${response.status}`);
        const data = await response.text();
        
        if (data.includes('ATOM') || data.includes('HETATM')) {
          return { data, format: 'pdb', source: 'Uploaded File' };
        }
        throw new Error('Invalid file format');
      } catch (err) {
        console.warn('File fetch failed:', err);
      }
    }

    // Raw PDB data
    if (trimmedInput.includes('ATOM') || trimmedInput.includes('HETATM')) {
      return { data: trimmedInput, format: 'pdb', source: 'Direct PDB Data' };
    }

    // PDB ID
    if (/^[1-9][A-Za-z0-9]{3}$/i.test(trimmedInput)) {
      try {
        console.log(`🧬 Fetching PDB ID: ${trimmedInput}`);
        const pdbUrl = `https://api.allorigins.win/raw?url=${encodeURIComponent(`https://files.rcsb.org/download/${trimmedInput.toUpperCase()}.pdb`)}`;
        const response = await fetch(pdbUrl);
        
        if (response.ok) {
          const data = await response.text();
          if (data.includes('ATOM') || data.includes('HETATM')) {
            return { data, format: 'pdb', source: 'RCSB PDB Database' };
          }
        }
      } catch (err) {
        console.warn('PDB fetch failed:', err);
      }
    }

    // PubChem
    try {
      console.log(`💊 Fetching from PubChem: ${trimmedInput}`);
      const isNumeric = /^[0-9]+$/.test(trimmedInput);
      const searchType = isNumeric ? 'cid' : 'name';
      
      const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/${searchType}/${encodeURIComponent(trimmedInput)}/SDF`;
      const response = await fetch(`https://api.allorigins.win/raw?url=${encodeURIComponent(pubchemUrl)}`);
      
      if (response.ok) {
        const data = await response.text();
        if (data.includes('M  END') && data.includes('V2000')) {
          return { data, format: 'sdf', source: 'PubChem Database' };
        }
      }
    } catch (err) {
      console.warn('PubChem fetch failed:', err);
    }

    throw new Error(`Could not resolve structure: "${trimmedInput}"`);
  };

  // Extract specific pose from multi-model PDB/PDBQT
  const extractPoseFromData = (data: string, poseIndex: number): string => {
    if (!isResultsMode || poseIndex === 0) return data;
    
    // For PDBQT files with multiple models
    if (data.includes('MODEL')) {
      const models = data.split(/MODEL\s+\d+/);
      if (models.length > poseIndex + 1) {
        const selectedModel = models[poseIndex + 1];
        const endIndex = selectedModel.indexOf('ENDMDL');
        if (endIndex !== -1) {
          return `MODEL     ${poseIndex + 1}\n${selectedModel.substring(0, endIndex)}\nENDMDL\nEND`;
        }
      }
    }
    
    return data; // Return original if pose extraction fails
  };

  // Add binding pockets visualization
  const addPocketVisualization = useCallback((viewer: any) => {
    if (!pocketsVisible || !viewer) return;

    console.log('🎯 Adding pocket visualization...');

    // Add primary pocket (red)
    if (pocketInfo) {
      const [x, y, z] = pocketInfo.center;
      const [sx, sy, sz] = pocketInfo.size;
      
      console.log(`Primary pocket: center(${x}, ${y}, ${z}), size(${sx}, ${sy}, ${sz})`);
      
      // Add wireframe box for binding pocket
      viewer.addBox({
        center: { x, y, z },
        dimensions: { w: sx, h: sy, d: sz },
        color: 'red',
        alpha: 0.4,
        wireframe: true
      });
      
      // Add center sphere
      viewer.addSphere({
        center: { x, y, z },
        radius: 3.0,
        color: 'red',
        alpha: 0.8
      });
    }

    // Add other pockets
    if (allPockets && allPockets.length > 1) {
      const colors = ['blue', 'green', 'orange', 'purple', 'yellow', 'cyan'];
      
      allPockets.slice(1).forEach((pocket, idx) => {
        const [x, y, z] = pocket.center;
        const [sx, sy, sz] = pocket.size;
        const color = colors[idx % colors.length];
        
        viewer.addBox({
          center: { x, y, z },
          dimensions: { w: sx, h: sy, d: sz },
          color: color,
          alpha: 0.3,
          wireframe: true
        });
        
        viewer.addSphere({
          center: { x, y, z },
          radius: 2.0,
          color: color,
          alpha: 0.6
        });
      });
    }
  }, [pocketInfo, allPockets, pocketsVisible]);

  // Main initialization effect
  useEffect(() => {
    if (!moleculePath || moleculePath.trim() === '') return;
    
    const initializeViewer = async () => {
      try {
        setIsLoading(true);
        setError(null);
        setHasViewer(false);
        setDataSource('');

        // Load 3Dmol script
        if (!window.$3Dmol) {
          await new Promise<void>((resolve, reject) => {
            const script = document.createElement('script');
            script.src = 'https://3Dmol.csb.pitt.edu/build/3Dmol-min.js';
            script.onload = () => {
              setTimeout(() => {
                if (window.$3Dmol) resolve();
                else reject(new Error('3Dmol library not available'));
              }, 200);
            };
            script.onerror = () => reject(new Error('Failed to load 3Dmol script'));
            document.head.appendChild(script);
          });
        }

        // Set up viewer element
        if (!viewerRef.current || !containerRef.current) {
          throw new Error('Viewer elements not ready');
        }

        const containerRect = containerRef.current.getBoundingClientRect();
        const minWidth = Math.max(containerRect.width, 300);
        const minHeight = Math.max(containerRect.height, 300);

        viewerRef.current.style.width = `${minWidth}px`;
        viewerRef.current.style.height = `${minHeight}px`;
        viewerRef.current.style.minWidth = '300px';
        viewerRef.current.style.minHeight = '300px';

        await new Promise(resolve => setTimeout(resolve, 100));

        const viewer = window.$3Dmol!.createViewer(viewerRef.current, {
          backgroundColor: 'white'
        });

        if (!viewer) {
          throw new Error('Viewer creation returned null');
        }

        viewer.setBackgroundColor('white');
        viewer.resize();
        viewerInstanceRef.current = viewer;

        console.log(`🔍 Resolving structure: "${moleculePath}" (Results mode: ${isResultsMode})`);
        const { data, format, source } = await resolveStructureData(moleculePath);
        setDataSource(source);
        setCurrentData(data);

        console.log(`✅ Successfully resolved from ${source}`);
        
        // Handle pose selection for docking results
        const poseData = extractPoseFromData(data, selectedPose);
        
        // Clear and load model
        viewer.clear();
        const model = viewer.addModel(poseData, format);
        if (!model) {
          throw new Error('Failed to create model from data');
        }

        // Apply styling
        if (moleculeType === 'complex' || isResultsMode) {
          console.log('🔬 Styling docked complex');
          
          // Protein: cartoon (blue-ish)
          viewer.setStyle({ chain: 'A' }, { 
            cartoon: { 
              color: 'lightblue', 
              alpha: 0.8 
            } 
          });
          
          // Ligand: ball-and-stick (colorful)
          viewer.setStyle({ hetflag: true }, { 
            stick: { 
              radius: 0.5, 
              colorscheme: 'Jmol' 
            },
            sphere: { 
              radius: 0.6, 
              colorscheme: 'Jmol' 
            }
          });
          
          // Highlight binding site
          if (pocketInfo) {
            const [x, y, z] = pocketInfo.center;
            const [sx, sy, sz] = pocketInfo.size;
            
            viewer.setStyle({
              and: [
                { atom: 'CA' },
                { x: { $gte: x - sx/2, $lte: x + sx/2 }},
                { y: { $gte: y - sy/2, $lte: y + sy/2 }},
                { z: { $gte: z - sz/2, $lte: z + sz/2 }}
              ]
            }, {
              cartoon: { color: 'orange', alpha: 0.9 }
            });
          }
          
        } else if (format === 'pdb' && (moleculeType === 'protein' || source.includes('PDB'))) {
          viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
        } else {
          viewer.setStyle({}, { 
            stick: { radius: 0.4, colorscheme: 'Jmol' },
            sphere: { radius: 0.4, colorscheme: 'Jmol' }
          });
        }

        // Add pocket visualizations
        addPocketVisualization(viewer);

        // ENHANCED ZOOM AND RESIZE
        await new Promise(resolve => setTimeout(resolve, 100));
        
        // First, zoom to all atoms
        viewer.zoomTo();
        
        // Wait a bit, then apply additional zoom for better visibility
        setTimeout(() => {
          if (moleculeType === 'complex' || isResultsMode) {
            // For docked complexes, zoom out a bit to see the full structure
            viewer.zoom(0.8, 1000);
          } else {
            // For individual molecules, zoom in for detail
            viewer.zoom(1.2, 1000);
          }
          viewer.render();
        }, 500);
        
        viewer.render();
        setHasViewer(true);
        console.log(`🎉 Successfully rendered: ${moleculePath} (Pose ${selectedPose + 1})`);

      } catch (err) {
        console.error('❌ Viewer initialization failed:', err);
        setError(err instanceof Error ? err.message : 'Failed to load viewer');
        setHasViewer(false);
        setDataSource('');
      } finally {
        setIsLoading(false);
      }
    };

    initializeViewer();
  }, [moleculePath, moleculeType, isResultsMode, addPocketVisualization]);

  // Handle pose changes
  useEffect(() => {
    if (!viewerInstanceRef.current || !currentData || !isResultsMode) return;
    
    console.log(`🔄 Switching to pose ${selectedPose + 1}`);
    
    const viewer = viewerInstanceRef.current;
    const poseData = extractPoseFromData(currentData, selectedPose);
    
    // Clear and reload with new pose
    viewer.clear();
    const model = viewer.addModel(poseData, 'pdbqt');
    
    if (model) {
      // Reapply styling
      viewer.setStyle({ chain: 'A' }, { 
        cartoon: { color: 'lightblue', alpha: 0.8 } 
      });
      
      viewer.setStyle({ hetflag: true }, { 
        stick: { radius: 0.5, colorscheme: 'Jmol' },
        sphere: { radius: 0.6, colorscheme: 'Jmol' }
      });
      
      // Re-add pockets
      addPocketVisualization(viewer);
      
      // Re-zoom and render
      viewer.zoomTo();
      setTimeout(() => {
        viewer.zoom(0.8, 500);
        viewer.render();
      }, 200);
    }
  }, [selectedPose, currentData, isResultsMode, addPocketVisualization]);

  // Handle pocket visibility toggle
  useEffect(() => {
    if (!viewerInstanceRef.current || !isResultsMode) return;
    
    const viewer = viewerInstanceRef.current;
    
    // Clear only shapes (pockets), keep molecules
    viewer.clear();
    
    // Re-add the current model
    const poseData = extractPoseFromData(currentData, selectedPose);
    const model = viewer.addModel(poseData, 'pdbqt');
    
    if (model) {
      // Reapply styling
      viewer.setStyle({ chain: 'A' }, { 
        cartoon: { color: 'lightblue', alpha: 0.8 } 
      });
      
      viewer.setStyle({ hetflag: true }, { 
        stick: { radius: 0.5, colorscheme: 'Jmol' },
        sphere: { radius: 0.6, colorscheme: 'Jmol' }
      });
      
      // Add pockets if visible
      if (pocketsVisible) {
        addPocketVisualization(viewer);
      }
      
      viewer.render();
    }
  }, [pocketsVisible, currentData, selectedPose, isResultsMode, addPocketVisualization]);

  // Toggle pocket visibility
  const togglePockets = () => {
    setPocketsVisible(!pocketsVisible);
  };

  // No molecule loaded state
  if (!moleculePath || moleculePath.trim() === '') {
    return (
      <div className="w-full h-full flex items-center justify-center" style={{ height }}>
        <div className="text-center text-gray-400">
          <div className="text-6xl mb-6">
            {moleculeType === 'protein' ? '🧬' : 
             moleculeType === 'ligand' ? '💊' : 
             moleculeType === 'complex' ? '🔬' : '⚗️'}
          </div>
          <p className="text-sm">No {moleculeType} loaded</p>
        </div>
      </div>
    );
  }

  const isFilePath = moleculePath.includes('/') || moleculePath.includes('\\');
  const displayName = isFilePath ? 
    (isResultsMode ? 'Docking Results' : 'Uploaded File') : 
    moleculePath;

  return (
    <div 
      ref={containerRef}
      className="w-full h-full relative bg-gradient-to-br from-gray-900 to-gray-700 rounded-lg overflow-hidden" 
      style={{ height, minHeight: '350px' }}
    >
      {/* 3D Viewer */}
      <div
        ref={viewerRef}
        className="absolute inset-0"
        style={{ 
          width: '100%', 
          height: '100%',
          minWidth: '300px',
          minHeight: '300px',
          visibility: hasViewer && !isLoading && !error ? 'visible' : 'hidden'
        }}
      />

      {/* Placeholder/Error state */}
      {(!hasViewer || isLoading || error) && (
        <div className="absolute inset-0">
          <div className="absolute inset-0 opacity-20">
            <div className="absolute top-4 left-4 w-2 h-2 bg-white rounded-full animate-pulse"></div>
            <div className="absolute top-8 right-8 w-1 h-1 bg-white rounded-full animate-pulse delay-300"></div>
            <div className="absolute bottom-6 left-8 w-1 h-1 bg-white rounded-full animate-pulse delay-500"></div>
            <div className="absolute bottom-4 right-4 w-2 h-2 bg-white rounded-full animate-pulse delay-700"></div>
            
            <svg className="w-full h-full opacity-10">
              <defs>
                <pattern id={`grid-${moleculeType}`} width="20" height="20" patternUnits="userSpaceOnUse">
                  <path d="M 20 0 L 0 0 0 20" fill="none" stroke="white" strokeWidth="0.5"/>
                </pattern>
              </defs>
              <rect width="100%" height="100%" fill={`url(#grid-${moleculeType})`} />
            </svg>
          </div>

          <div className="absolute inset-0 flex items-center justify-center">
            <div className="text-center p-6">
              <div className="relative mb-4">
                <div 
                  className="text-6xl transform transition-transform duration-1000 hover:scale-110" 
                  style={{ color }}
                >
                  {moleculeType === 'protein' ? '🧬' : 
                   moleculeType === 'ligand' ? '💊' : 
                   moleculeType === 'complex' ? '🔬' : '⚗️'}
                </div>
                
                <div className="absolute inset-0 flex items-center justify-center">
                  <div className="w-16 h-16 border border-white opacity-30 rounded-full animate-spin"></div>
                </div>
                <div className="absolute inset-0 flex items-center justify-center">
                  <div className="w-20 h-20 border border-white opacity-20 rounded-full animate-ping"></div>
                </div>
              </div>

              <div className="bg-black bg-opacity-50 rounded-lg p-4 backdrop-blur-sm">
                <div className="text-lg font-semibold text-white mb-2">
                  {displayName}
                </div>
                <div className="text-sm text-gray-300 mb-1 capitalize">
                  {moleculeType === 'complex' ? 'Docked Complex' : `${moleculeType} Structure`}
                  {selectedPose > 0 && ` - Pose ${selectedPose + 1}`}
                </div>
                <div className="text-xs text-gray-400 font-mono bg-gray-800 px-2 py-1 rounded max-w-full truncate">
                  {error ? `Error: ${error}` : 
                   isLoading ? `Loading ${isResultsMode ? 'docking results' : 'structure'}...` :
                   dataSource ? `Source: ${dataSource}` :
                   moleculePath}
                </div>
                
                {isResultsMode && pocketInfo && (
                  <div className="mt-2 text-xs text-gray-300">
                    <div>Binding Site: ({pocketInfo.center.map(c => c.toFixed(1)).join(', ')})</div>
                    <div>Method: {pocketInfo.method} • Confidence: {pocketInfo.confidence}</div>
                  </div>
                )}
                
                <div className="mt-4 flex justify-center space-x-2">
                  <div className={`w-2 h-2 rounded-full ${
                    error ? 'bg-red-400' : 
                    isLoading ? 'bg-yellow-400 animate-pulse' : 
                    hasViewer ? 'bg-green-400' : 'bg-gray-400'
                  }`}></div>
                  <div className="w-2 h-2 bg-blue-400 rounded-full animate-pulse delay-200"></div>
                  <div className="w-2 h-2 bg-purple-400 rounded-full animate-pulse delay-400"></div>
                </div>
                <div className="text-xs text-gray-400 mt-2">
                  {error ? 'Visualization failed' : 
                   isLoading ? 'Loading molecular data...' : 
                   hasViewer ? (isResultsMode ? 'Docking complex rendered' : `Rendered from ${dataSource}`) : 'Initializing...'}
                </div>
              </div>

              <div className="flex justify-center space-x-2 mt-4">
                <button 
                  className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors"
                  onClick={() => window.location.reload()}
                >
                  {error ? 'Retry' : 'Rotate'}
                </button>
                <button className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors">
                  Zoom
                </button>
                <button className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors">
                  Reset
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Status indicator */}
      <div className="absolute bottom-4 right-4 bg-black bg-opacity-70 rounded px-3 py-2">
        <div className="flex items-center space-x-2">
          <div className={`w-3 h-3 rounded-full ${
            error ? 'bg-red-400' : 
            isLoading ? 'bg-yellow-400 animate-pulse' : 
            hasViewer ? 'bg-green-400' : 'bg-gray-400'
          }`}></div>
          <span className="text-xs text-white">
            {error ? 'Error' : 
             isLoading ? 'Loading' : 
             hasViewer ? (isResultsMode ? 'Docked' : '3D Active') : 'Rendering'}
          </span>
        </div>
      </div>

      {/* Pocket toggle for results mode */}
      {isResultsMode && hasViewer && (pocketInfo || allPockets) && (
        <div className="absolute top-4 right-4 bg-black bg-opacity-70 rounded px-3 py-2">
          <div className="flex items-center space-x-2">
            <span className="text-xs text-white">Pockets:</span>
            <button
              onClick={togglePockets}
              className={`text-xs px-2 py-1 rounded transition-colors ${
                pocketsVisible ? 'bg-red-500 text-white hover:bg-red-600' : 'bg-gray-500 text-white hover:bg-gray-600'
              }`}
            >
              {pocketsVisible ? 'Hide' : 'Show'}
            </button>
          </div>
        </div>
      )}

      {/* Pose info overlay */}
      {isResultsMode && hasViewer && poseData && poseData[selectedPose] && (
        <div className="absolute top-4 left-4 bg-black bg-opacity-70 rounded px-3 py-2">
          <div className="text-xs text-white">
            <div>Pose {selectedPose + 1}/{poseData.length}</div>
            <div>Affinity: {poseData[selectedPose].affinity} kcal/mol</div>
          </div>
        </div>
      )}
    </div>
  );
};
