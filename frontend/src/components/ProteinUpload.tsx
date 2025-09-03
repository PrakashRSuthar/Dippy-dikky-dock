// src/components/ProteinUpload.tsx
import { useState } from 'react';
import { Settings, Activity, Upload, Atom } from 'lucide-react';
import { useApiHealth } from '../hooks/useApiHealth';
import { useFileUpload } from '../hooks/useFileUpload';
import { MoleculeVisualization } from './MoleculeVisualization';

interface ProteinUploadProps {
  onMoleculesPaired: (proteinInput: string, ligandInput: string) => void;
}

export const ProteinUpload = ({ onMoleculesPaired }: ProteinUploadProps) => {
  const [showSettings, setShowSettings] = useState(false);
  
  // Protein state
  const [activeProteinTab, setActiveProteinTab] = useState<'database' | 'upload'>('database');
  const [proteinId, setProteinId] = useState('');
  const [proteinInput, setProteinInput] = useState<string>(''); // This goes to backend
  const [proteinSource, setProteinSource] = useState<'database' | 'upload' | null>(null);
  
  // Ligand state
  const [activeLigandTab, setActiveLigandTab] = useState<'database' | 'upload'>('database');
  const [ligandId, setLigandId] = useState('');
  const [ligandInput, setLigandInput] = useState<string>(''); // This goes to backend
  const [ligandSource, setLigandSource] = useState<'database' | 'upload' | null>(null);

  const { isHealthy, isChecking } = useApiHealth();
  const { uploadFile, isUploading, uploadedFileName } = useFileUpload();

  // Protein handlers
  const handleProteinFileUpload = async (file: File) => {
    const result = await uploadFile(file, 'protein');
    if (result?.file_path) {
      setProteinInput(result.file_path); // Backend expects file_path
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
      setProteinInput(proteinId.trim()); // Backend expects ID string
      setProteinSource('database');
    }
  };

  // Ligand handlers  
  const handleLigandFileUpload = async (file: File) => {
    const result = await uploadFile(file, 'ligand');
    if (result?.file_path) {
      setLigandInput(result.file_path); // Backend expects file_path
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
      setLigandInput(ligandId.trim()); // Backend expects ID string
      setLigandSource('database');
    }
  };

  // Continue to docking when both molecules are ready
  const canProceed = proteinInput && ligandInput;
  const handleProceedToDocking = () => {
    if (canProceed) {
      onMoleculesPaired(proteinInput, ligandInput);
    }
  };

  const isValidProteinId = proteinId.trim().length >= 4;
  const isValidLigandId = ligandId.trim().length >= 2;

  return (
    <div className="max-w-7xl mx-auto p-6">
      {/* Header */}
      <div className="flex items-center justify-between mb-8">
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
              {isChecking ? 'Checking...' : isHealthy ? 'System OK' : 'System Error'}
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

      {/* Steps */}
      <div className="flex items-center space-x-8 mb-8">
        <div className="flex items-center space-x-2">
          <div className="w-8 h-8 bg-blue-600 rounded-full flex items-center justify-center">
            <Upload className="w-4 h-4 text-white" />
          </div>
          <span className="font-medium text-blue-600">Upload</span>
        </div>
        <div className="flex-1 h-px bg-gray-300" />
        <div className="flex items-center space-x-2">
          <div className="w-8 h-8 bg-gray-300 rounded-full flex items-center justify-center">
            <span className="text-xs">⚡</span>
          </div>
          <span className="text-gray-500">Docking</span>
        </div>
        <div className="flex-1 h-px bg-gray-300" />
        <div className="flex items-center space-x-2">
          <div className="w-8 h-8 bg-gray-300 rounded-full flex items-center justify-center">
            <span className="text-xs">📊</span>
          </div>
          <span className="text-gray-500">Results</span>
        </div>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Left Column: Protein & Ligand Upload */}
        <div className="space-y-6">
          {/* Protein Section */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center space-x-2">
                <span>🧬</span>
                <span>Protein Structure</span>
                {proteinInput && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-600 mt-1">Upload protein structure or fetch from PDB/AlphaFold</p>
            </div>

            <div className="flex border-b">
              <button
                onClick={() => setActiveProteinTab('database')}
                className={`px-6 py-3 text-sm font-medium border-b-2 ${
                  activeProteinTab === 'database'
                    ? 'text-blue-600 border-blue-600 bg-blue-50'
                    : 'text-gray-600 border-transparent hover:text-gray-900'
                }`}
              >
                Database
              </button>
              <button
                onClick={() => setActiveProteinTab('upload')}
                className={`px-6 py-3 text-sm font-medium border-b-2 ${
                  activeProteinTab === 'upload'
                    ? 'text-blue-600 border-blue-600 bg-blue-50'
                    : 'text-gray-600 border-transparent hover:text-gray-900'
                }`}
              >
                Upload
              </button>
            </div>

            <div className="p-6">
              {activeProteinTab === 'database' && (
                <div>
                  <input
                    type="text"
                    value={proteinId}
                    onChange={(e) => setProteinId(e.target.value.toUpperCase())}
                    placeholder="1HSG or P12345"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md text-center font-mono focus:outline-none focus:ring-2 focus:ring-blue-500"
                  />
                  <div className="text-xs text-gray-500 mt-2 text-center">
                    PDB ID (1HSG) or UniProt ID (P12345)
                  </div>
                  <button
                    onClick={handleProteinDatabase}
                    disabled={!isValidProteinId}
                    className={`w-full mt-4 py-2 px-4 rounded-md font-medium transition-colors ${
                      isValidProteinId
                        ? 'bg-blue-600 text-white hover:bg-blue-700'
                        : 'bg-gray-300 text-gray-500 cursor-not-allowed'
                    }`}
                  >
                    Fetch Protein
                  </button>
                </div>
              )}

              {activeProteinTab === 'upload' && (
                <div
                  className="border-2 border-dashed border-gray-300 rounded-xl p-8 text-center hover:border-gray-400 transition-colors"
                  onDrop={handleProteinDrop}
                  onDragOver={(e) => e.preventDefault()}
                >
                  <Upload className="w-8 h-8 text-gray-400 mx-auto mb-2" />
                  <p className="text-sm text-gray-600 mb-2">Drop PDB/PDBQT files here</p>
                  <input
                    type="file"
                    accept=".pdb,.pdbqt"
                    onChange={handleProteinFileSelect}
                    className="hidden"
                    id="protein-file"
                  />
                  <label
                    htmlFor="protein-file"
                    className="inline-flex items-center px-3 py-1 border border-gray-300 rounded-md bg-white text-sm cursor-pointer hover:bg-gray-50"
                  >
                    Browse
                  </label>
                  {isUploading && <div className="mt-2 text-xs text-blue-600">Uploading...</div>}
                </div>
              )}

              {proteinInput && (
                <div className="mt-4 p-3 bg-green-50 rounded-lg">
                  <div className="text-sm text-green-800 flex items-center space-x-2">
                    <span>✅</span>
                    <span>
                      Protein ready: {proteinSource === 'upload' ? uploadedFileName : proteinInput}
                    </span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Ligand Section */}
          <div className="bg-white rounded-xl border shadow-sm">
            <div className="px-6 py-4 border-b">
              <h2 className="text-lg font-semibold flex items-center space-x-2">
                <span>💊</span>
                <span>Ligand Structure</span>
                {ligandInput && <span className="text-green-600">✅</span>}
              </h2>
              <p className="text-sm text-gray-600 mt-1">Upload ligand or fetch from PubChem database</p>
            </div>

            <div className="flex border-b">
              <button
                onClick={() => setActiveLigandTab('database')}
                className={`px-6 py-3 text-sm font-medium border-b-2 ${
                  activeLigandTab === 'database'
                    ? 'text-green-600 border-green-600 bg-green-50'
                    : 'text-gray-600 border-transparent hover:text-gray-900'
                }`}
              >
                PubChem
              </button>
              <button
                onClick={() => setActiveLigandTab('upload')}
                className={`px-6 py-3 text-sm font-medium border-b-2 ${
                  activeLigandTab === 'upload'
                    ? 'text-green-600 border-green-600 bg-green-50'
                    : 'text-gray-600 border-transparent hover:text-gray-900'
                }`}
              >
                Upload
              </button>
            </div>

            <div className="p-6">
              {activeLigandTab === 'database' && (
                <div>
                  <input
                    type="text"
                    value={ligandId}
                    onChange={(e) => setLigandId(e.target.value)}
                    placeholder="Aspirin or 2244"
                    className="w-full px-3 py-2 border border-gray-300 rounded-md text-center font-mono focus:outline-none focus:ring-2 focus:ring-green-500"
                  />
                  <div className="text-xs text-gray-500 mt-2 text-center">
                    Compound name (Aspirin) or PubChem CID (2244)
                  </div>
                  <button
                    onClick={handleLigandDatabase}
                    disabled={!isValidLigandId}
                    className={`w-full mt-4 py-2 px-4 rounded-md font-medium transition-colors ${
                      isValidLigandId
                        ? 'bg-green-600 text-white hover:bg-green-700'
                        : 'bg-gray-300 text-gray-500 cursor-not-allowed'
                    }`}
                  >
                    Fetch Ligand
                  </button>
                </div>
              )}

              {activeLigandTab === 'upload' && (
                <div
                  className="border-2 border-dashed border-gray-300 rounded-xl p-8 text-center hover:border-gray-400 transition-colors"
                  onDrop={handleLigandDrop}
                  onDragOver={(e) => e.preventDefault()}
                >
                  <Atom className="w-8 h-8 text-gray-400 mx-auto mb-2" />
                  <p className="text-sm text-gray-600 mb-2">Drop SDF/MOL/MOL2 files here</p>
                  <input
                    type="file"
                    accept=".sdf,.mol,.mol2"
                    onChange={handleLigandFileSelect}
                    className="hidden"
                    id="ligand-file"
                  />
                  <label
                    htmlFor="ligand-file"
                    className="inline-flex items-center px-3 py-1 border border-gray-300 rounded-md bg-white text-sm cursor-pointer hover:bg-gray-50"
                  >
                    Browse
                  </label>
                </div>
              )}

              {ligandInput && (
                <div className="mt-4 p-3 bg-green-50 rounded-lg">
                  <div className="text-sm text-green-800 flex items-center space-x-2">
                    <span>✅</span>
                    <span>
                      Ligand ready: {ligandSource === 'upload' ? uploadedFileName : ligandInput}
                    </span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Proceed Button */}
          <button
            onClick={handleProceedToDocking}
            disabled={!canProceed}
            className={`w-full py-4 px-6 rounded-xl font-semibold text-lg transition-colors ${
              canProceed
                ? 'bg-blue-600 text-white hover:bg-blue-700 shadow-lg'
                : 'bg-gray-300 text-gray-500 cursor-not-allowed'
            }`}
          >
            {canProceed ? '⚡ Start Molecular Docking' : 'Select Protein & Ligand to Continue'}
          </button>
        </div>

        {/* Right Column: 3D Visualization */}
        <div className="bg-white rounded-xl border shadow-sm">
          <div className="px-6 py-4 border-b">
            <h2 className="text-lg font-semibold flex items-center space-x-2">
              <span>🔬</span>
              <span>3D Structure Preview</span>
            </h2>
            <p className="text-sm text-gray-600 mt-1">Interactive molecular visualization</p>
          </div>
          <div className="p-6">
            <div className="h-96 bg-gray-50 rounded-lg flex items-center justify-center">
              {proteinInput || ligandInput ? (
                <MoleculeVisualization 
                  proteinPath={proteinInput}
                  ligandPath={ligandInput}
                />
              ) : (
                <div className="text-center text-gray-500">
                  <div className="w-16 h-16 mx-auto mb-4 opacity-30">🧬</div>
                  <p>Select a protein or ligand to view 3D structure</p>
                </div>
              )}
            </div>
          </div>
        </div>
      </div>

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
