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
      setError(null);
      
      const apiBase = import.meta.env.VITE_API_BASE || 'http://localhost:8000';
      
      // Fetch results
      const resultsResponse = await fetch(`${apiBase}/api/jobs/${jobId}/result`);
      if (!resultsResponse.ok) {
        const errorData = await resultsResponse.json();
        throw new Error(errorData.detail || 'Failed to fetch results');
      }
      const resultsData = await resultsResponse.json();
      setResults(resultsData);
      
      // Fetch file list
      try {
        const filesResponse = await fetch(`${apiBase}/api/jobs/${jobId}/files`);
        if (filesResponse.ok) {
          const filesData = await filesResponse.json();
          setFiles(filesData.files || []);
        }
      } catch (err) {
        console.warn('Failed to fetch file list:', err);
      }
      
      // Simulate pose data (you can replace this with actual CSV parsing)
      const mockPoses = Array.from({ length: resultsData.total_poses }, (_, i) => ({
        pose: i + 1,
        affinity: (resultsData.best_affinity + Math.random() * 2 * (i / resultsData.total_poses)).toFixed(2),
        rmsd_lb: (Math.random() * 2).toFixed(2),
        rmsd_ub: (Math.random() * 3 + 1).toFixed(2),
      }));
      setPoseData(mockPoses);
      
    } catch (err) {
      setError(err instanceof Error ? err.message : 'An unknown error occurred');
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
          <p className="text-gray-600">Loading docking results for Job ID: {jobId}</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center p-4">
        <div className="text-center bg-white p-8 rounded-lg shadow-md border border-red-200">
          <XCircle className="w-16 h-16 text-red-500 mx-auto mb-4" />
          <h2 className="text-xl font-bold text-gray-800 mb-2">Error Loading Results</h2>
          <p className="text-gray-600 mb-6 max-w-md">{error}</p>
          <div className="flex gap-4 justify-center">
            <button
              onClick={onBack}
              className="px-4 py-2 bg-gray-600 text-white rounded-md hover:bg-gray-700"
            >
              Go Back
            </button>
            <button
              onClick={fetchResults}
              className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700"
            >
              Retry
            </button>
          </div>
        </div>
      </div>
    );
  }

  if (!results) {
    return (
      <div className="min-h-screen bg-gray-50 flex items-center justify-center">
        <div className="text-center">
          <p className="text-gray-600">No results available for this job.</p>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-50">
      {/* Header */}
      <div className="bg-white border-b sticky top-0 z-10">
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
                <p className="text-sm text-gray-600 font-mono">Job ID: {jobId}</p>
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
        <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-4 gap-4 mb-6">
          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Best Affinity</p>
                <p className="text-2xl font-bold text-green-600">{results.best_affinity}</p>
                <p className="text-xs text-gray-500">kcal/mol</p>
              </div>
              <Zap className="w-8 h-8 text-green-600 opacity-70" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Total Poses</p>
                <p className="text-2xl font-bold text-blue-600">{results.total_poses}</p>
                <p className="text-xs text-gray-500">conformations found</p>
              </div>
              <Atom className="w-8 h-8 text-blue-600 opacity-70" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Binding Sites</p>
                <p className="text-2xl font-bold text-purple-600">{results.all_pockets?.length || 1}</p>
                <p className="text-xs text-gray-500">detected</p>
              </div>
              <BarChart3 className="w-8 h-8 text-purple-600 opacity-70" />
            </div>
          </div>

          <div className="bg-white rounded-xl border shadow-sm p-4">
            <div className="flex items-center justify-between">
              <div>
                <p className="text-sm text-gray-600">Primary Pocket Size</p>
                <p className="text-lg font-bold text-orange-600">
                  {results.pocket_size[0].toFixed(1)}×{results.pocket_size[1].toFixed(1)}×{results.pocket_size[2].toFixed(1)}
                </p>
                <p className="text-xs text-gray-500">Ångströms</p>
              </div>
              <div className="text-3xl opacity-70">📦</div>
            </div>
          </div>
        </div>
        
        <div className="space-y-6">
          {/* TOP SECTION: Visualization + Pockets */}
          <div className="grid grid-cols-1 xl:grid-cols-3 gap-6">
            
            {/* 3D Visualization */}
            <div className="xl:col-span-2 bg-white rounded-xl border shadow-sm">
              <div className="px-6 py-4 border-b flex justify-between items-center">
                <h2 className="text-lg font-semibold flex items-center space-x-2">
                  <span>🔬</span>
                  <span>Docked Complex</span>
                  <span className="text-sm bg-green-100 text-green-800 px-2 py-1 rounded-full">
                    Pose {selectedPose + 1}
                  </span>
                </h2>
              </div>
              <div className="p-4 md:p-6">
                <div className="h-[50vh] min-h-[500px] border border-gray-200 rounded-lg bg-gray-50">
                  <MoleculeVisualization 
                    // ✅ FIX: Switched to proteinPath and dockedPath for results mode
                    proteinPath={results.files?.prepared_protein}
                    dockedPath={results.files?.docked_poses}
                    height="100%"
                    pocketInfo={{
                      center: results.pocket_center,
                      size: results.pocket_size,
                    }}
                    allPockets={results.all_pockets}
                    selectedPose={selectedPose}
                    showPockets={true}
                    isResultsMode={true}
                    poseData={poseData}
                  />
                </div>
              </div>
            </div>

            {/* Binding Pockets List */}
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
              <div className="p-4 flex-1 overflow-auto" style={{ maxHeight: 'calc(50vh + 100px)' }}>
                <div className="space-y-3">
                  {results.all_pockets?.map((pocket, idx) => (
                    <div key={idx} className={`p-3 border rounded-lg transition-colors hover:bg-gray-50 ${idx === 0 ? 'bg-blue-50 border-blue-200' : 'border-gray-200'}`}>
                      <div className="flex items-center justify-between mb-2">
                        <div className="font-medium text-sm">
                          {idx === 0 && <span className="text-blue-600">★ </span>}
                          Pocket {idx + 1}
                        </div>
                        <span className={`px-2 py-1 rounded-full text-xs font-medium ${
                          pocket.confidence === 'high' ? 'bg-green-100 text-green-800' :
                          pocket.confidence === 'medium' ? 'bg-yellow-100 text-yellow-800' :
                          'bg-gray-100 text-gray-800'
                        }`}>
                          {pocket.confidence}
                        </span>
                      </div>
                      <div className="text-xs text-gray-600 space-y-1">
                        <div><strong>Center:</strong> ({pocket.center?.map((c: number) => c.toFixed(1)).join(', ')})</div>
                        <div><strong>Size:</strong> {pocket.size?.map((s: number) => s.toFixed(1)).join('×')} Å</div>
                        <div><strong>Method:</strong> {pocket.method}</div>
                      </div>
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>

          {/* BOTTOM SECTION: Chart and Table */}
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
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
        </div>
      </div>
    </div>
  );
};