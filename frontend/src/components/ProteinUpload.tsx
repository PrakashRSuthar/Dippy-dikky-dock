// src/components/ProteinUpload.tsx
import { useState } from 'react';
import { Settings, Activity, Upload, Atom } from 'lucide-react';
import { useApiHealth } from '../hooks/useApiHealth';
import { useFileUpload } from '../hooks/useFileUpload';
import { useDocking } from '../hooks/useDocking';
import { MoleculeVisualization } from './MoleculeVisualization';

interface ProteinUploadProps {
  onDockingStarted: (jobId: string) => void;
}

export const ProteinUpload = ({ onDockingStarted }: ProteinUploadProps) => {
  const [showSettings, setShowSettings] = useState(false);
  
  // Protein state
  const [activeProteinTab, setActiveProteinTab] = useState<'database' | 'upload'>('database');
  const [proteinId, setProteinId] = useState('1HSG'); // Pre-fill example
  const [proteinInput, setProteinInput] = useState<string>('');
  const [proteinSource, setProteinSource] = useState<'database' | 'upload' | null>(null);
  
  // Ligand state
  const [activeLigandTab, setActiveLigandTab] = useState<'database' | 'upload'>('database');
  const [ligandId, setLigandId] = useState('aspirin'); // Pre-fill example
  const [ligandInput, setLigandInput] = useState<string>('');
  const [ligandSource, setLigandSource] = useState<'database' | 'upload' | null>(null);
  
  // Retention policy state
  const [retentionPolicy, setRetentionPolicy] = useState<'save' | 'delete' | 'keep7d'>('keep7d');

  const { isHealthy, isChecking } = useApiHealth();
  const { uploadFile, isUploading, uploadedFileName } = useFileUpload();
  const { startDocking, isStarting } = useDocking();

  // Protein handlers
  const handleProteinFileUpload = async (file: File) => {
    if (!isHealthy) {
      alert('Backend not connected. Please start your FastAPI server.');
      return;
    }
    const result = await uploadFile(file, 'protein');
    if (result?.file_path) {
      setProteinInput(result.file_path);
      setProteinSource('upload');
    }
  };

  const handleProteinFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) handleProteinFileUpload(file);
  };

  const handleProteinDrop = (e: React.DragEvent) => {
    e.preventDefault();
    const file = e.dataTransfer.files[0];
    if (file && (file.name.endsWith('.pdb') || file.name.endsWith('.pdbqt'))) {
      handleProteinFileUpload(file);
    }
  };

  const handleProteinDatabase = () => {
    if (proteinId.trim().length >= 4) {
      setProteinInput(proteinId.trim());
      setProteinSource('database');
    }
  };

  // Ligand handlers  
  const handleLigandFileUpload = async (file: File) => {
    if (!isHealthy) {
      alert('Backend not connected. Please start your FastAPI server.');
      return;
    }
    const result = await uploadFile(file, 'ligand');
    if (result?.file_path) {
      setLigandInput(result.file_path);
      setLigandSource('upload');
    }
  };

  const handleLigandFileSelect = (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (file) handleLigandFileUpload(file);
  };

  const handleLigandDrop = (e: React.DragEvent) => {
    e.preventDefault();
    const file = e.dataTransfer.files[0];
    if (file && (file.name.endsWith('.sdf') || file.name.endsWith('.mol') || file.name.endsWith('.mol2'))) {
      handleLigandFileUpload(file);
    }
  };

  const handleLigandDatabase = () => {
    if (ligandId.trim().length >= 2) {
      setLigandInput(ligandId.trim());
      setLigandSource('database');
    }
  };

  // Start Docking Process
  const handleStartDocking = async () => {
    if (!proteinInput || !ligandInput) return;
    
    const dockingRequest = {
      protein_input: proteinInput,
      ligand_input: ligandInput,
      job_name: `${proteinInput}_${ligandInput}_${Date.now()}`,
      cleaning_policy: {
        keep_waters: false,
        keep_ions: true,
        keep_solvents: false,
        keep_cofactors: true
      },
      retention: retentionPolicy  // Include retention policy
    };

    const jobId = await startDocking(dockingRequest);
    if (jobId) {
      onDockingStarted(jobId);
    }
  };

  const canProceed = proteinInput && ligandInput && isHealthy;
  const isValidProteinId = proteinId.trim().length >= 4;
  const isValidLigandId = ligandId.trim().length >= 2;

  return (
    <div className="max-w-full mx-auto p-4 lg:p-6">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div className="flex items-center space-x-3">
          <div className="w-10 h-10 bg-blue-600 rounded-lg flex items-center justify-center">
            <Activity className="w-6 h-6 text-white" />
          </div>
          <div>
            <h1 className="text-2xl font-bold">MolecularDock Pro</h1>
            <p className="text-sm text-gray-600">Professional Molecular Docking Platform</p>
          </div>
          <div className="flex items-center space-x-2 ml-4">
            <div className={`w-3 h-3 rounded-full ${
              isChecking ? 'bg-yellow-400 animate-pulse' : isHealthy ? 'bg-green-400' : 'bg-red-400'
            }`} />
            <span className="text-xs text-gray-600">
              {isChecking ? 'Checking...' : isHealthy ? 'System OK' : 'Backend Offline'}
            </span>
          </div>
        </div>
        <button
          onClick={() => setShowSettings(true)}
          className="flex items-center space-x-2 px-4 py-2 text-gray-600 hover:bg-gray-100 rounded-lg"
        >
          <Settings className="w-4 h-4" />
          <span>Settings</span>
        </button>
      </div>

      {/* Backend Status */}
      {!isHealthy && (
        <div className="mb-6 p-4 bg-yellow-50 border border-yellow-200 rounded-lg">
          <div className="flex items-center space-x-2">
            <span className="text-yellow-600">ℹ️</span>
            <div>
              <h3 className="text-yellow-800 font-medium">Backend Required for Docking</h3>
              <p className="text-yellow-700 text-sm">
                Start backend: <code className="bg-yellow-100 px-1 rounded">python api/docking_api.py</code>
              </p>
            </div>
          </div>
        </div>
      )}

      {/* Main Content - Full Width Layout */}
      <div className="grid grid-cols-1 xl:grid-cols-2 gap-6">
        
        {/* Left Column: Protein Section */}
        <div className="space-y-6">
          {/* Protein Upload - More Compact */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center space-x-2">
                <span>🧬</span>
                <span>Protein Structure</span>
                {proteinInput && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-600 mt-1">Upload or fetch from database</p>
            </div>

            <div className="flex border-b">
              <button
                onClick={() => setActiveProteinTab('database')}
                className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 ${
                  activeProteinTab === 'database'
                    ? 'text-blue-600 border-blue-600 bg-blue-50'
                    : 'text-gray-600 border-transparent'
                }`}
              >
                Database
              </button>
              <button
                onClick={() => setActiveProteinTab('upload')}
                className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 ${
                  activeProteinTab === 'upload'
                    ? 'text-blue-600 border-blue-600 bg-blue-50'
                    : 'text-gray-600 border-transparent'
                }`}
              >
                Upload
              </button>
            </div>

            <div className="p-4">
              {activeProteinTab === 'database' && (
                <div className="space-y-3">
                  <input
                    type="text"
                    value={proteinId}
                    onChange={(e) => setProteinId(e.target.value.toUpperCase())}
                    placeholder="1HSG or P12345"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md text-center font-mono focus:outline-none focus:ring-2 focus:ring-blue-500"
                  />
                  <div className="text-xs text-gray-500 text-center">PDB ID or UniProt ID</div>
                  <button
                    onClick={handleProteinDatabase}
                    disabled={!isValidProteinId}
                    className={`w-full py-2 px-4 rounded-md font-medium transition-colors ${
                      isValidProteinId
                        ? 'bg-blue-600 text-white hover:bg-blue-700'
                        : 'bg-gray-300 text-gray-500 cursor-not-allowed'
                    }`}
                  >
                    Use Protein ID
                  </button>
                </div>
              )}

              {activeProteinTab === 'upload' && (
                <div>
                  <div
                    className={`border-2 border-dashed rounded-xl p-6 text-center transition-colors ${
                      isHealthy 
                        ? 'border-gray-300 hover:border-gray-400' 
                        : 'border-gray-200 bg-gray-50'
                    }`}
                    onDrop={handleProteinDrop}
                    onDragOver={(e) => e.preventDefault()}
                  >
                    <Upload className={`w-6 h-6 mx-auto mb-2 ${isHealthy ? 'text-gray-400' : 'text-gray-300'}`} />
                    <p className="text-sm text-gray-600 mb-2">Drop PDB/PDBQT files here</p>
                    <input
                      type="file"
                      accept=".pdb,.pdbqt"
                      onChange={handleProteinFileSelect}
                      className="hidden"
                      id="protein-file"
                      disabled={!isHealthy}
                    />
                    <label
                      htmlFor="protein-file"
                      className={`inline-flex items-center px-3 py-1 border border-gray-300 rounded-md bg-white text-sm ${
                        isHealthy ? 'cursor-pointer hover:bg-gray-50' : 'cursor-not-allowed opacity-50'
                      }`}
                    >
                      Browse Files
                    </label>
                    {isUploading && <div className="mt-2 text-xs text-blue-600">Uploading...</div>}
                    {!isHealthy && <div className="mt-2 text-xs text-red-600">Backend required for uploads</div>}
                  </div>
                </div>
              )}

              {proteinInput && (
                <div className="mt-4 p-3 bg-green-50 rounded-lg">
                  <div className="text-sm text-green-800 flex items-center space-x-2">
                    <span>✅</span>
                    <span>Ready: {proteinSource === 'upload' ? uploadedFileName : proteinInput}</span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Protein Visualization - EXPANDED */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h3 className="text-lg font-semibold flex items-center space-x-2">
                <span>🔬</span>
                <span>Protein Preview</span>
              </h3>
              <p className="text-sm text-gray-600">3D structure visualization</p>
            </div>
            <div className="p-4">
              {/* EXPANDED HEIGHT - Better fit for protein structures */}
              <div className="h-96 md:h-[500px] lg:h-[600px] border border-gray-200 rounded-lg bg-gray-50">
                <MoleculeVisualization 
                  moleculePath={proteinInput}
                  moleculeType="protein"
                  color="#3b82f6"
                  height="100%"
                />
              </div>
              
              {/* Protein Info Panel */}
              {proteinInput && (
                <div className="mt-4 grid grid-cols-1 sm:grid-cols-2 gap-4 text-sm">
                  <div className="bg-blue-50 p-3 rounded-lg">
                    <div className="font-medium text-blue-800">Source</div>
                    <div className="text-blue-700">
                      {proteinSource === 'database' ? 'RCSB PDB Database' : 'Uploaded File'}
                    </div>
                  </div>
                  <div className="bg-blue-50 p-3 rounded-lg">
                    <div className="font-medium text-blue-800">Structure ID</div>
                    <div className="text-blue-700 font-mono">
                      {proteinSource === 'upload' ? uploadedFileName : proteinInput}
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>

        {/* Right Column: Ligand Section */}
        <div className="space-y-6">
          {/* Ligand Upload - More Compact */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center space-x-2">
                <span>💊</span>
                <span>Ligand Structure</span>
                {ligandInput && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-600 mt-1">Upload or fetch from PubChem</p>
            </div>

            <div className="flex border-b">
              <button
                onClick={() => setActiveLigandTab('database')}
                className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 ${
                  activeLigandTab === 'database'
                    ? 'text-green-600 border-green-600 bg-green-50'
                    : 'text-gray-600 border-transparent'
                }`}
              >
                PubChem
              </button>
              <button
                onClick={() => setActiveLigandTab('upload')}
                className={`flex-1 px-4 py-3 text-sm font-medium border-b-2 ${
                  activeLigandTab === 'upload'
                    ? 'text-green-600 border-green-600 bg-green-50'
                    : 'text-gray-600 border-transparent'
                }`}
              >
                Upload
              </button>
            </div>

            <div className="p-4">
              {activeLigandTab === 'database' && (
                <div className="space-y-3">
                  <input
                    type="text"
                    value={ligandId}
                    onChange={(e) => setLigandId(e.target.value)}
                    placeholder="aspirin or 2244"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md text-center font-mono focus:outline-none focus:ring-2 focus:ring-green-500"
                  />
                  <div className="text-xs text-gray-500 text-center">Compound name or CID</div>
                  <button
                    onClick={handleLigandDatabase}
                    disabled={!isValidLigandId}
                    className={`w-full py-2 px-4 rounded-md font-medium transition-colors ${
                      isValidLigandId
                        ? 'bg-green-600 text-white hover:bg-green-700'
                        : 'bg-gray-300 text-gray-500 cursor-not-allowed'
                    }`}
                  >
                    Use Ligand ID
                  </button>
                </div>
              )}

              {activeLigandTab === 'upload' && (
                <div>
                  <div
                    className={`border-2 border-dashed rounded-xl p-6 text-center transition-colors ${
                      isHealthy 
                        ? 'border-gray-300 hover:border-gray-400' 
                        : 'border-gray-200 bg-gray-50'
                    }`}
                    onDrop={handleLigandDrop}
                    onDragOver={(e) => e.preventDefault()}
                  >
                    <Atom className={`w-6 h-6 mx-auto mb-2 ${isHealthy ? 'text-gray-400' : 'text-gray-300'}`} />
                    <p className="text-sm text-gray-600 mb-2">Drop SDF/MOL/MOL2 files here</p>
                    <input
                      type="file"
                      accept=".sdf,.mol,.mol2"
                      onChange={handleLigandFileSelect}
                      className="hidden"
                      id="ligand-file"
                      disabled={!isHealthy}
                    />
                    <label
                      htmlFor="ligand-file"
                      className={`inline-flex items-center px-3 py-1 border border-gray-300 rounded-md bg-white text-sm ${
                        isHealthy ? 'cursor-pointer hover:bg-gray-50' : 'cursor-not-allowed opacity-50'
                      }`}
                    >
                      Browse Files
                    </label>
                    {!isHealthy && <div className="mt-2 text-xs text-red-600">Backend required for uploads</div>}
                  </div>
                </div>
              )}

              {ligandInput && (
                <div className="mt-4 p-3 bg-green-50 rounded-lg">
                  <div className="text-sm text-green-800 flex items-center space-x-2">
                    <span>✅</span>
                    <span>Ready: {ligandSource === 'upload' ? uploadedFileName : ligandInput}</span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Ligand Visualization - EXPANDED */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h3 className="text-lg font-semibold flex items-center space-x-2">
                <span>⚗️</span>
                <span>Ligand Preview</span>
              </h3>
              <p className="text-sm text-gray-600">Chemical structure visualization</p>
            </div>
            <div className="p-4">
              {/* EXPANDED HEIGHT - Better fit for ligand structures */}
              <div className="h-96 md:h-[500px] lg:h-[600px] border border-gray-200 rounded-lg bg-gray-50">
                <MoleculeVisualization 
                  moleculePath={ligandInput}
                  moleculeType="ligand"
                  color="#10b981"
                  height="100%"
                />
              </div>

              {/* Ligand Info Panel */}
              {ligandInput && (
                <div className="mt-4 grid grid-cols-1 sm:grid-cols-2 gap-4 text-sm">
                  <div className="bg-green-50 p-3 rounded-lg">
                    <div className="font-medium text-green-800">Source</div>
                    <div className="text-green-700">
                      {ligandSource === 'database' ? 'PubChem Database' : 'Uploaded File'}
                    </div>
                  </div>
                  <div className="bg-green-50 p-3 rounded-lg">
                    <div className="font-medium text-green-800">Compound</div>
                    <div className="text-green-700 font-mono">
                      {ligandSource === 'upload' ? uploadedFileName : ligandInput}
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>
        </div>
      </div>

      {/* Data Retention Policy Selection */}
      {canProceed && (
        <div className="mt-8 bg-white rounded-xl border shadow-sm p-6">
          <h3 className="text-lg font-semibold mb-4 flex items-center space-x-2">
            <span>💾</span>
            <span>Data Retention Policy</span>
          </h3>
          <p className="text-sm text-gray-600 mb-4">
            Choose what to do with generated files after docking completes:
          </p>
          
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <label className="flex items-start space-x-3 cursor-pointer p-4 border-2 rounded-lg transition-colors hover:bg-gray-50">
              <input
                type="radio"
                name="retention"
                value="save"
                checked={retentionPolicy === 'save'}
                onChange={(e) => setRetentionPolicy(e.target.value as any)}
                className="mt-1"
              />
              <div>
                <div className="font-medium flex items-center space-x-2">
                  <span>💾</span>
                  <span>Save Permanently</span>
                </div>
                <div className="text-sm text-gray-600 mt-1">
                  Move all files to permanent storage in data/results/
                </div>
              </div>
            </label>
            
            <label className="flex items-start space-x-3 cursor-pointer p-4 border-2 rounded-lg transition-colors hover:bg-gray-50">
              <input
                type="radio"
                name="retention"
                value="keep7d"
                checked={retentionPolicy === 'keep7d'}
                onChange={(e) => setRetentionPolicy(e.target.value as any)}
                className="mt-1"
              />
              <div>
                <div className="font-medium flex items-center space-x-2">
                  <span>⏳</span>
                  <span>Keep Temporarily</span>
                </div>
                <div className="text-sm text-gray-600 mt-1">
                  Keep files in temp directory, auto-delete after 7 days
                </div>
              </div>
            </label>
            
            <label className="flex items-start space-x-3 cursor-pointer p-4 border-2 rounded-lg transition-colors hover:bg-gray-50">
              <input
                type="radio"
                name="retention"
                value="delete"
                checked={retentionPolicy === 'delete'}
                onChange={(e) => setRetentionPolicy(e.target.value as any)}
                className="mt-1"
              />
              <div>
                <div className="font-medium flex items-center space-x-2">
                  <span>🗑️</span>
                  <span>Delete Immediately</span>
                </div>
                <div className="text-sm text-gray-600 mt-1">
                  Remove all generated files after results are processed
                </div>
              </div>
            </label>
          </div>
        </div>
      )}

      {/* Start Docking Button */}
      {canProceed && (
        <div className="mt-8 flex justify-center">
          <button
            onClick={handleStartDocking}
            disabled={!canProceed || isStarting}
            className={`px-8 py-4 rounded-lg font-semibold text-lg transition-colors flex items-center space-x-3 ${
              canProceed && !isStarting
                ? 'bg-blue-600 text-white hover:bg-blue-700 shadow-lg'
                : 'bg-gray-300 text-gray-500 cursor-not-allowed'
            }`}
          >
            {isStarting ? (
              <>
                <div className="w-5 h-5 border-2 border-white border-t-transparent rounded-full animate-spin"></div>
                <span>Starting Docking...</span>
              </>
            ) : (
              <>
                <span>⚡</span>
                <span>Start Molecular Docking</span>
              </>
            )}
          </button>
        </div>
      )}

      {/* Settings Modal */}
      {showSettings && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
          <div className="bg-white rounded-lg p-6 w-96">
            <h3 className="text-lg font-semibold mb-4">Settings</h3>
            <div className="mb-4">
              <label className="block text-sm font-medium mb-1">API Base URL</label>
              <input
                type="url"
                defaultValue={import.meta.env.VITE_API_BASE || 'http://localhost:8000'}
                className="w-full px-3 py-2 border rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
              />
            </div>
            <div className="flex space-x-3">
              <button
                onClick={() => setShowSettings(false)}
                className="flex-1 px-4 py-2 border rounded-md hover:bg-gray-50"
              >
                Cancel
              </button>
              <button
                onClick={() => setShowSettings(false)}
                className="flex-1 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
              >
                Save
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};
