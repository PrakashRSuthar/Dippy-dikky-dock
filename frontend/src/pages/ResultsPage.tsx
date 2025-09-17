// src/pages/ResultsPage.tsx
import { useState, useEffect } from 'react';
import { ArrowLeft, Download, RefreshCw, BarChart3, Zap, Atom, CheckCircle, XCircle } from 'lucide-react';
import { MoleculeVisualization } from '../components/MoleculeVisualization';
import { ResultsTable } from '../components/ResultsTable';
import { ResultsChart } from '../components/ResultsChart';

interface ResultsPageProps {
  jobId: string;
  onBack: () => void;
}

interface PocketInfo {
  center: number[];
  size: number[];
  confidence: string;
  method: string;
}

interface FileInfo {
  name: string;
  filename: string;
  path: string;
  size: string;
  download_url: string;
}

interface DockingResults {
  job_id: string;
  status: string;
  best_affinity: number;
  total_poses: number;
  pocket_center: number[];
  pocket_size: number[];
  method: string;
  confidence: string;
  run_directory: string;
  all_pockets: PocketInfo[];
  files: Record<string, string>;
  protein_analysis: string;
  cleaning_policy: any;
  retention: string;
  permanent_location?: string;
}

export const ResultsPage = ({ jobId, onBack }: ResultsPageProps) => {
  const [results, setResults] = useState<DockingResults | null>(null);
  const [files, setFiles] = useState<FileInfo[]>([]);
  const [poseData, setPoseData] = useState<any[]>([]);
  const [selectedPose, setSelectedPose] = useState<number>(0);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetchResults();
  }, [jobId]);

  const fetchResults = async () => {
    try {
      setLoading(true);
      
      // Fetch results
      const resultsResponse = await fetch(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/jobs/${jobId}/result`);
      
      if (!resultsResponse.ok) {
        throw new Error('Failed to fetch results');
      }

      const resultsData = await resultsResponse.json();
      setResults(resultsData);
      
      // Fetch file list
      try {
        const filesResponse = await fetch(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/jobs/${jobId}/files`);
        if (filesResponse.ok) {
          const filesData = await filesResponse.json();
          setFiles(filesData.files || []);
        }
      } catch (err) {
        console.warn('Failed to fetch file list:', err);
      }
      
      // Simulate pose data (in real implementation, this would come from your CSV)
      const mockPoses = Array.from({ length: resultsData.total_poses }, (_, i) => ({
        pose: i + 1,
        affinity: (resultsData.best_affinity + Math.random() * 3).toFixed(2),
        rmsd_lb: (Math.random() * 2).toFixed(2),
        rmsd_ub: (Math.random() * 3 + 1).toFixed(2),
      }));
      setPoseData(mockPoses);
      
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load results');
    } finally {
      setLoading(false);
    }
  };

  const downloadFile = (url: string, filename: string) => {
    const link = document.createElement('a');
    link.href = `${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}${url}`;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const downloadAllFiles = () => {
    files.forEach((file, index) => {
      setTimeout(() => {
        downloadFile(file.download_url, file.filename);
      }, index * 200); // Stagger downloads
    });
  };

  if (loading) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="text-center">
          <RefreshCw className="w-8 h-8 animate-spin text-blue-600 mx-auto mb-4" />
          <p className="text-gray-600">Loading docking results...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="text-center">
          <XCircle className="w-16 h-16 text-red-600 mx-auto mb-4" />
          <h2 className="text-xl font-bold text-gray-800 mb-2">Error Loading Results</h2>
          <p className="text-gray-600 mb-4">{error}</p>
          <button
            onClick={onBack}
            className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
          >
            Go Back
          </button>
        </div>
      </div>
    );
  }

  if (!results) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="text-center">
          <p className="text-gray-600">No results available</p>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white border-b">
        <div className="max-w-full mx-auto px-6 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-4">
              <button
                onClick={onBack}
                className="flex items-center space-x-2 text-gray-600 hover:text-gray-900 transition-colors"
              >
                <ArrowLeft className="w-5 h-5" />
                <span>Back</span>
              </button>
              <div>
                <h1 className="text-2xl font-bold text-gray-900 flex items-center space-x-2">
                  <CheckCircle className="w-6 h-6 text-green-600" />
                  <span>Docking Results</span>
                </h1>
                <p className="text-sm text-gray-600">Job ID: {jobId}</p>
              </div>
            </div>
            <div className="flex items-center space-x-3">
              {files.length > 0 && (
                <button
                  onClick={downloadAllFiles}
                  className="flex items-center space-x-2 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
                >
                  <Download className="w-4 h-4" />
                  <span>Download All</span>
                </button>
              )}
            </div>
          </div>
        </div>
      </div>

      <div className="max-w-full mx-auto px-4 py-6">
        {/* Summary Cards */}
        <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mb-6">
          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Best Affinity</p>
                <p className="text-2xl font-bold text-green-600">{results.best_affinity}</p>
                <p className="text-xs text-gray-500">kcal/mol</p>
              </div>
              <Zap className="w-8 h-8 text-green-600" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Total Poses</p>
                <p className="text-2xl font-bold text-blue-600">{results.total_poses}</p>
                <p className="text-xs text-gray-500">conformations</p>
              </div>
              <Atom className="w-8 h-8 text-blue-600" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Binding Sites</p>
                <p className="text-2xl font-bold text-purple-600">{results.all_pockets?.length || 1}</p>
                <p className="text-xs text-gray-500">detected</p>
              </div>
              <BarChart3 className="w-8 h-8 text-purple-600" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Primary Pocket</p>
                <p className="text-lg font-bold text-orange-600">
                  {results.pocket_size[0].toFixed(1)}×{results.pocket_size[1].toFixed(1)}×{results.pocket_size[2].toFixed(1)}
                </p>
                <p className="text-xs text-gray-500">Ångströms³</p>
              </div>
              <div className="w-8 h-8 text-orange-600 text-2xl">📦</div>
            </div>
          </div>
        </div>

        {/* ENHANCED LAYOUT: Full Width Top Section + Two Column Bottom */}
        <div className="space-y-6">
          
          {/* TOP SECTION: Large 3D Visualization + Binding Pockets Table */}
          <div className="grid grid-cols-1 xl:grid-cols-3 gap-6">
            
            {/* EXPANDED 3D Visualization - Takes 2/3 width on XL screens */}
            <div className="xl:col-span-2 bg-white rounded-xl border shadow-sm">
              <div className="px-6 py-4 border-b">
                <h2 className="text-lg font-semibold flex items-center space-x-2">
                  <span>🔬</span>
                  <span>Docked Complex</span>
                  <span className="text-sm bg-green-100 text-green-800 px-2 py-1 rounded-full">
                    Pose {selectedPose + 1}
                  </span>
                </h2>
              </div>
              <div className="p-6">
                {/* MUCH LARGER VISUALIZATION - Increased from h-80 to h-[600px] */}
                <div className="h-[600px] border border-gray-200 rounded-lg bg-gray-50">
                  <MoleculeVisualization 
                    moleculePath={results.files?.docked_poses}
                    moleculeType="complex"
                    color="#8b5cf6"
                    height="100%"
                    pocketInfo={{
                      center: results.pocket_center,
                      size: results.pocket_size,
                      confidence: results.confidence,
                      method: results.method
                    }}
                    allPockets={results.all_pockets}
                    selectedPose={selectedPose}
                    showPockets={true}
                    jobId={jobId}
                    isResultsMode={true}
                    poseData={poseData}
                  />
                </div>
                
                {/* Enhanced Pose Controls */}
                <div className="mt-6 grid grid-cols-1 md:grid-cols-3 gap-4">
                  {/* Pose Selector */}
                  <div className="flex items-center space-x-2">
                    <span className="text-sm text-gray-600 font-medium">Pose:</span>
                    <select
                      value={selectedPose}
                      onChange={(e) => setSelectedPose(Number(e.target.value))}
                      className="flex-1 border border-gray-300 rounded-md px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
                    >
                      {poseData.map((_, idx) => (
                        <option key={idx} value={idx}>
                          Pose {idx + 1} ({poseData[idx]?.affinity} kcal/mol)
                        </option>
                      ))}
                    </select>
                  </div>
                  
                  {/* Navigation Buttons */}
                  <div className="flex items-center space-x-2">
                    <button
                      onClick={() => setSelectedPose(Math.max(0, selectedPose - 1))}
                      disabled={selectedPose === 0}
                      className="flex-1 px-4 py-2 border border-gray-300 rounded-md text-sm hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                    >
                      ← Previous
                    </button>
                    <button
                      onClick={() => setSelectedPose(Math.min(poseData.length - 1, selectedPose + 1))}
                      disabled={selectedPose === poseData.length - 1}
                      className="flex-1 px-4 py-2 border border-gray-300 rounded-md text-sm hover:bg-gray-50 disabled:opacity-50 disabled:cursor-not-allowed transition-colors"
                    >
                      Next →
                    </button>
                  </div>

                  {/* Current Pose Info */}
                  <div className="text-center">
                    <div className="text-sm text-gray-600">Current Affinity</div>
                    <div className="text-lg font-bold text-green-600">
                      {poseData[selectedPose]?.affinity} kcal/mol
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Binding Pockets Table - Takes 1/3 width */}
            <div className="bg-white rounded-xl border shadow-sm flex flex-col">
              <div className="px-6 py-4 border-b">
                <h2 className="text-lg font-semibold flex items-center space-x-2">
                  <span>🎯</span>
                  <span>Binding Pockets</span>
                  <span className="text-sm bg-purple-100 text-purple-800 px-2 py-1 rounded-full">
                    {results.all_pockets?.length || 0}
                  </span>
                </h2>
              </div>
              <div className="p-4 flex-1 overflow-auto">
                <div className="space-y-3">
                  {results.all_pockets?.map((pocket, idx) => (
                    <div key={idx} className={`p-3 border rounded-lg transition-colors hover:bg-gray-50 ${
                      idx === 0 ? 'bg-blue-50 border-blue-200' : 'border-gray-200'
                    }`}>
                      <div className="flex items-center justify-between mb-2">
                        <div className="font-medium text-sm">
                          {idx === 0 && <span className="text-blue-600">★ </span>}
                          Pocket {idx + 1}
                        </div>
                        <span className={`px-2 py-1 rounded-full text-xs ${
                          pocket.confidence === 'high' ? 'bg-green-100 text-green-800' :
                          pocket.confidence === 'medium' ? 'bg-yellow-100 text-yellow-800' :
                          'bg-gray-100 text-gray-800'
                        }`}>
                          {pocket.confidence}
                        </span>
                      </div>
                      <div className="text-xs text-gray-600 space-y-1">
                        <div><strong>Center:</strong> ({pocket.center?.map((c: number) => c.toFixed(1)).join(', ')})</div>
                        <div><strong>Size:</strong> {pocket.size?.map((s: number) => s.toFixed(1)).join('×')} Ų</div>
                        <div><strong>Method:</strong> {pocket.method}</div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>

          {/* BOTTOM SECTION: Two Column Layout */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            
            {/* Left Column */}
            <div className="space-y-6">
              
              {/* Binding Affinity Chart */}
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="text-lg font-semibold flex items-center space-x-2">
                    <span>📊</span>
                    <span>Binding Affinity Distribution</span>
                  </h2>
                </div>
                <div className="p-6">
                  <div className="h-64">
                    <ResultsChart data={poseData} selectedPose={selectedPose} />
                  </div>
                </div>
              </div>

              {/* Results Table */}
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="text-lg font-semibold flex items-center space-x-2">
                    <span>📋</span>
                    <span>Docking Poses</span>
                  </h2>
                </div>
                <div className="p-6">
                  <ResultsTable 
                    data={poseData} 
                    selectedPose={selectedPose}
                    onPoseSelect={setSelectedPose}
                  />
                </div>
              </div>
            </div>

            {/* Right Column */}
            <div className="space-y-6">
              
              {/* Generated Files */}
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="text-lg font-semibold flex items-center space-x-2">
                    <span>📁</span>
                    <span>Generated Files</span>
                    <span className="text-sm bg-green-100 text-green-800 px-2 py-1 rounded-full">
                      {files.length} files
                    </span>
                  </h2>
                </div>
                <div className="p-6">
                  {files.length > 0 ? (
                    <div className="space-y-3">
                      {files.map((file, idx) => (
                        <div key={idx} className="flex items-center justify-between p-3 border border-gray-200 rounded-lg hover:bg-gray-50 transition-colors">
                          <div className="flex items-center space-x-3 min-w-0 flex-1">
                            <div className="w-10 h-10 bg-blue-100 rounded-lg flex items-center justify-center flex-shrink-0">
                              <span className="text-sm">📄</span>
                            </div>
                            <div className="min-w-0 flex-1">
                              <div className="font-medium text-sm">{file.name}</div>
                              <div className="text-xs text-gray-500">{file.filename} • {file.size}</div>
                              <div className="text-xs text-gray-400 font-mono truncate">
                                {file.path}
                              </div>
                            </div>
                          </div>
                          <button
                            onClick={() => downloadFile(file.download_url, file.filename)}
                            className="px-3 py-1 bg-blue-600 text-white text-xs rounded-md hover:bg-blue-700 transition-colors flex-shrink-0 ml-3"
                          >
                            Download
                          </button>
                        </div>
                      ))}
                    </div>
                  ) : (
                    <div className="text-center py-8 text-gray-500">
                      <div className="text-4xl mb-2">📄</div>
                      <p>No files available for download</p>
                    </div>
                  )}
                  
                  {files.length > 0 && (
                    <div className="mt-4 pt-4 border-t">
                      <button
                        onClick={downloadAllFiles}
                        className="w-full px-4 py-2 bg-gray-100 text-gray-700 rounded-md hover:bg-gray-200 transition-colors"
                      >
                        📦 Download All Files
                      </button>
                    </div>
                  )}
                </div>
              </div>

              {/* Data Retention Status */}
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="text-lg font-semibold flex items-center space-x-2">
                    <span>💾</span>
                    <span>Data Retention</span>
                  </h2>
                </div>
                <div className="p-6">
                  <div className={`p-4 rounded-lg border-2 ${
                    results.retention === 'save' ? 'bg-green-50 border-green-200' :
                    results.retention === 'delete' ? 'bg-red-50 border-red-200' :
                    'bg-yellow-50 border-yellow-200'
                  }`}>
                    <div className="flex items-center space-x-2">
                      <span className="text-2xl">
                        {results.retention === 'save' ? '💾' : results.retention === 'delete' ? '🗑️' : '⏳'}
                      </span>
                      <div>
                        <div className="font-medium">
                          {results.retention === 'save' ? 'Saved Permanently' :
                           results.retention === 'delete' ? 'Will Be Deleted' :
                           'Temporary Storage (7 days)'}
                        </div>
                        <div className="text-sm text-gray-600">
                          {results.retention === 'save' ? 'Files are stored permanently in data/results/' :
                           results.retention === 'delete' ? 'All files will be deleted after processing' :
                           'Files will be automatically deleted in 7 days'}
                        </div>
                      </div>
                    </div>
                  </div>
                  
                  {results.permanent_location && (
                    <div className="mt-3 p-3 bg-blue-50 rounded-lg">
                      <div className="text-sm">
                        <strong>Permanent Location:</strong>
                        <div className="font-mono text-xs mt-1 text-blue-800 break-all">
                          {results.permanent_location}
                        </div>
                      </div>
                    </div>
                  )}
                </div>
              </div>

              {/* Analysis Summary */}
              <div className="bg-white rounded-xl border shadow-sm">
                <div className="px-6 py-4 border-b">
                  <h2 className="text-lg font-semibold flex items-center space-x-2">
                    <span>🧪</span>
                    <span>Analysis Summary</span>
                  </h2>
                </div>
                <div className="p-6 space-y-4">
                  <div>
                    <h3 className="font-medium text-gray-900 mb-2">Protein Analysis</h3>
                    <p className="text-sm text-gray-600">{results.protein_analysis}</p>
                  </div>
                  
                  <div>
                    <h3 className="font-medium text-gray-900 mb-2">Primary Binding Site</h3>
                    <div className="text-sm text-gray-600 space-y-1">
                      <p><strong>Center:</strong> ({results.pocket_center.map(c => c.toFixed(1)).join(', ')})</p>
                      <p><strong>Size:</strong> {results.pocket_size.map(s => s.toFixed(1)).join(' × ')} Ų</p>
                      <p><strong>Method:</strong> {results.method}</p>
                      <p><strong>Confidence:</strong> {results.confidence}</p>
                    </div>
                  </div>

                  <div>
                    <h3 className="font-medium text-gray-900 mb-2">Cleaning Policy</h3>
                    <div className="text-sm text-gray-600">
                      <div className="grid grid-cols-2 gap-2">
                        <span>Waters: {results.cleaning_policy?.keep_waters ? '✅ Kept' : '❌ Removed'}</span>
                        <span>Ions: {results.cleaning_policy?.keep_ions ? '✅ Kept' : '❌ Removed'}</span>
                        <span>Solvents: {results.cleaning_policy?.keep_solvents ? '✅ Kept' : '❌ Removed'}</span>
                        <span>Cofactors: {results.cleaning_policy?.keep_cofactors ? '✅ Kept' : '❌ Removed'}</span>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};
