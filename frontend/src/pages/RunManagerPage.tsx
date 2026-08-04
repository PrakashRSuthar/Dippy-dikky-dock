import { useState, useEffect, useCallback, useMemo } from 'react';
import {
  Search, RefreshCw, Trash2, Eye, FolderOpen, ChevronRight,
  ChevronDown, File, BarChart3, Calendar, HardDrive, XCircle,
  Folder, X, Trophy, Database
} from 'lucide-react';

const apiBase = ((import.meta as unknown as Record<string, Record<string, string>>).env ?? {}).VITE_API_BASE || 'http://localhost:8000';

interface RunFile {
  name: string;
  size: number;
  modified: string;
}

interface RunInfo {
  name: string;
  date: string;
  file_count: number;
  total_size: number;
  has_result: boolean;
  files: RunFile[];
}

interface BenchmarkResult {
  name: string;
  date: string;
  file_count: number;
  total_size: number;
}

type RunStatus = 'completed' | 'running' | 'failed' | 'unknown';

function formatFileSize(bytes: number): string {
  if (bytes === 0) return '0 B';
  const k = 1024;
  const sizes = ['B', 'KB', 'MB', 'GB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  return `${parseFloat((bytes / Math.pow(k, i)).toFixed(1))} ${sizes[i]}`;
}

function formatDate(dateStr: string): string {
  const date = new Date(dateStr);
  return date.toLocaleDateString('en-US', {
    month: 'short',
    day: 'numeric',
    year: 'numeric',
    hour: '2-digit',
    minute: '2-digit'
  });
}

function getStatusFromRun(run: RunInfo): RunStatus {
  if (run.has_result) return 'completed';
  return 'unknown';
}

const statusConfig: Record<RunStatus, { color: string; bg: string; label: string }> = {
  completed: { color: 'text-green-400', bg: 'bg-green-900/30 border-green-700', label: 'Completed' },
  running: { color: 'text-yellow-400', bg: 'bg-yellow-900/30 border-yellow-700', label: 'Running' },
  failed: { color: 'text-red-400', bg: 'bg-red-900/30 border-red-700', label: 'Failed' },
  unknown: { color: 'text-gray-400', bg: 'bg-gray-800 border-gray-600', label: 'Pending' },
};

function StatusBadge({ status }: { status: RunStatus }) {
  const config = statusConfig[status];
  return (
    <span className={`inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium border ${config.bg} ${config.color}`}>
      {config.label}
    </span>
  );
}

const RunManagerPage = () => {
  const [runs, setRuns] = useState<RunInfo[]>([]);
  const [permanentRuns, setPermanentRuns] = useState<RunInfo[]>([]);
  const [benchmarks, setBenchmarks] = useState<BenchmarkResult[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');
  const [expandedRun, setExpandedRun] = useState<string | null>(null);
  const [runFiles, setRunFiles] = useState<Record<string, RunFile[]>>({});
  const [loadingFiles, setLoadingFiles] = useState<string | null>(null);
  const [deleteConfirmRun, setDeleteConfirmRun] = useState<string | null>(null);
  const [deleting, setDeleting] = useState<string | null>(null);
  const [toast, setToast] = useState<{ message: string; type: 'success' | 'error' } | null>(null);

  const showToast = useCallback((message: string, type: 'success' | 'error' = 'success') => {
    setToast({ message, type });
    setTimeout(() => setToast(null), 3000);
  }, []);

  const fetchRuns = useCallback(async () => {
    try {
      setLoading(true);
      setError(null);
      const res = await fetch(`${apiBase}/api/runs`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      setRuns(Array.isArray(data.temporary_runs) ? data.temporary_runs : (Array.isArray(data.runs) ? data.runs : []));
      setPermanentRuns(Array.isArray(data.permanent_runs) ? data.permanent_runs : []);

      // Fetch benchmarks separately
      try {
        const benchRes = await fetch(`${apiBase}/api/runs/benchmark_results`);
        if (benchRes.ok) {
          const benchData = await benchRes.json();
          setBenchmarks(benchData.files || benchData || []);
        }
      } catch {
        setBenchmarks([]);
      }
    } catch (e: unknown) {
      setError(e instanceof Error ? e.message : 'Failed to load runs');
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchRuns();
  }, [fetchRuns]);

  const toggleExpand = useCallback(async (runName: string) => {
    if (expandedRun === runName) {
      setExpandedRun(null);
      return;
    }
    setExpandedRun(runName);
    if (!runFiles[runName]) {
      setLoadingFiles(runName);
      try {
        const res = await fetch(`${apiBase}/api/runs/${encodeURIComponent(runName)}`);
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        const data = await res.json();
        setRunFiles(prev => ({ ...prev, [runName]: data.files || [] }));
      } catch {
        setRunFiles(prev => ({ ...prev, [runName]: [] }));
      } finally {
        setLoadingFiles(null);
      }
    }
  }, [expandedRun, runFiles]);

  const handleDelete = useCallback(async (runName: string) => {
    setDeleting(runName);
    try {
      const res = await fetch(`${apiBase}/api/runs/${encodeURIComponent(runName)}`, { method: 'DELETE' });
      if (!res.ok) throw new Error('Failed to delete');
      setRuns(prev => prev.filter(r => r.name !== runName));
      setPermanentRuns(prev => prev.filter(r => r.name !== runName));
      setDeleteConfirmRun(null);
      setExpandedRun(null);
      showToast('Run deleted successfully');
    } catch {
      showToast('Failed to delete run', 'error');
    } finally {
      setDeleting(null);
    }
  }, [showToast]);

  const handleViewResults = useCallback((runName: string) => {
    window.location.href = `${apiBase}/api/runs/${encodeURIComponent(runName)}/result.json`;
  }, []);

  const filteredRuns = useMemo(() => {
    if (!searchQuery) return runs;
    const q = searchQuery.toLowerCase();
    return runs.filter(r => r.name.toLowerCase().includes(q));
  }, [runs, searchQuery]);

  const filteredPermanentRuns = useMemo(() => {
    if (!searchQuery) return permanentRuns;
    const q = searchQuery.toLowerCase();
    return permanentRuns.filter(r => r.name.toLowerCase().includes(q));
  }, [permanentRuns, searchQuery]);

  const renderRunCard = (run: RunInfo) => {
    const status = getStatusFromRun(run);
    const isExpanded = expandedRun === run.name;
    const files = runFiles[run.name];

    return (
      <div
        key={run.name}
        className="bg-gray-800 rounded-xl border border-gray-700 overflow-hidden"
      >
        {/* Run card header */}
        <button
          onClick={() => toggleExpand(run.name)}
          className="w-full flex items-center gap-4 p-4 hover:bg-gray-750 transition-colors text-left"
        >
          <div className="flex-shrink-0 text-gray-400">
            {isExpanded ? (
              <ChevronDown className="w-5 h-5" />
            ) : (
              <ChevronRight className="w-5 h-5" />
            )}
          </div>
          <Folder className="w-5 h-5 text-blue-400 flex-shrink-0" />
          <div className="flex-1 min-w-0">
            <div className="flex items-center gap-3 flex-wrap">
              <span className="font-medium text-white truncate">{run.name}</span>
              <StatusBadge status={status} />
            </div>
            <div className="flex items-center gap-4 mt-1 text-xs text-gray-400 flex-wrap">
              <span className="flex items-center gap-1">
                <Calendar className="w-3 h-3" />
                {formatDate(run.date)}
              </span>
              <span className="flex items-center gap-1">
                <File className="w-3 h-3" />
                {run.file_count} files
              </span>
              <span className="flex items-center gap-1">
                <HardDrive className="w-3 h-3" />
                {formatFileSize(run.total_size)}
              </span>
            </div>
          </div>
          <div className="flex items-center gap-1 flex-shrink-0" onClick={(e) => e.stopPropagation()}>
            {run.has_result && (
              <button
                onClick={() => handleViewResults(run.name)}
                className="p-2 text-gray-400 hover:text-green-400 hover:bg-gray-700 rounded transition-colors"
                title="View Results"
              >
                <Eye className="w-4 h-4" />
              </button>
            )}
            {deleteConfirmRun === run.name ? (
              <div className="flex items-center gap-1">
                <button
                  onClick={() => handleDelete(run.name)}
                  disabled={deleting === run.name}
                  className="px-2 py-1 bg-red-600 text-white text-xs rounded hover:bg-red-700 disabled:opacity-50"
                >
                  {deleting === run.name ? '...' : 'Yes'}
                </button>
                <button
                  onClick={() => setDeleteConfirmRun(null)}
                  className="px-2 py-1 bg-gray-600 text-white text-xs rounded hover:bg-gray-500"
                >
                  No
                </button>
              </div>
            ) : (
              <button
                onClick={() => setDeleteConfirmRun(run.name)}
                className="p-2 text-gray-400 hover:text-red-400 hover:bg-gray-700 rounded transition-colors"
                title="Delete"
              >
                <Trash2 className="w-4 h-4" />
              </button>
            )}
          </div>
        </button>

        {/* Expanded file list */}
        {isExpanded && (
          <div className="border-t border-gray-700 bg-gray-900/50">
            {loadingFiles === run.name ? (
              <div className="flex items-center justify-center py-6">
                <RefreshCw className="w-5 h-5 animate-spin text-blue-400" />
              </div>
            ) : files && files.length > 0 ? (
              <div className="max-h-64 overflow-y-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-gray-700">
                      <th className="px-4 py-2 text-left text-xs font-semibold text-gray-400">Filename</th>
                      <th className="px-4 py-2 text-left text-xs font-semibold text-gray-400">Size</th>
                      <th className="px-4 py-2 text-left text-xs font-semibold text-gray-400">Modified</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-gray-800">
                    {files.map((f) => (
                      <tr key={f.name} className="hover:bg-gray-800/50">
                        <td className="px-4 py-2 text-gray-300 font-mono text-xs">{f.name}</td>
                        <td className="px-4 py-2 text-gray-400 text-xs">{formatFileSize(f.size)}</td>
                        <td className="px-4 py-2 text-gray-400 text-xs">{formatDate(f.modified)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            ) : (
              <div className="flex items-center justify-center py-6 text-gray-500 text-sm">
                No files in this run
              </div>
            )}
          </div>
        )}
      </div>
    );
  };

  return (
    <div className="min-h-screen bg-gray-900 text-gray-100">
      {toast && (
        <div className={`fixed top-4 right-4 z-50 px-4 py-3 rounded-lg shadow-lg transition-all ${
          toast.type === 'success' ? 'bg-green-600 text-white' : 'bg-red-600 text-white'
        }`}>
          {toast.message}
        </div>
      )}

      {/* Header */}
      <div className="bg-gray-800 border-b border-gray-700 sticky top-0 z-10">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <FolderOpen className="w-6 h-6 text-blue-400" />
              <h1 className="text-2xl font-bold">Run Manager</h1>
            </div>
            <button
              onClick={fetchRuns}
              disabled={loading}
              className="flex items-center gap-2 px-4 py-2 bg-gray-700 hover:bg-gray-600 rounded-lg transition-colors disabled:opacity-50"
            >
              <RefreshCw className={`w-4 h-4 ${loading ? 'animate-spin' : ''}`} />
              Refresh
            </button>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6 space-y-6">
        {/* Benchmark Section */}
        {benchmarks.length > 0 && (
          <div className="bg-gray-800 rounded-xl border border-gray-700 overflow-hidden">
            <div className="px-4 py-3 bg-gray-750 border-b border-gray-700 flex items-center gap-2">
              <Trophy className="w-4 h-4 text-yellow-400" />
              <h2 className="text-sm font-semibold text-gray-200 uppercase tracking-wider">Benchmark Results</h2>
            </div>
            <div className="p-4 grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
              {benchmarks.map((bench) => (
                <div
                  key={bench.name}
                  className="flex items-center justify-between p-3 bg-gray-900 rounded-lg border border-gray-700 hover:border-gray-600 transition-colors"
                >
                  <div className="flex items-center gap-3 min-w-0">
                    <BarChart3 className="w-4 h-4 text-yellow-400 flex-shrink-0" />
                    <div className="min-w-0">
                      <p className="text-sm font-medium text-white truncate">{bench.name}</p>
                      <p className="text-xs text-gray-400">{bench.file_count} files</p>
                    </div>
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}

        {/* Search */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-4">
          <div className="relative">
            <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 w-4 h-4 text-gray-400" />
            <input
              type="text"
              placeholder="Search runs..."
              value={searchQuery}
              onChange={(e) => setSearchQuery(e.target.value)}
              className="w-full pl-10 pr-4 py-2 bg-gray-700 border border-gray-600 rounded-lg text-white placeholder-gray-400 focus:outline-none focus:ring-2 focus:ring-blue-500 focus:border-transparent"
            />
            {searchQuery && (
              <button
                onClick={() => setSearchQuery('')}
                className="absolute right-3 top-1/2 transform -translate-y-1/2 text-gray-400 hover:text-white"
              >
                <X className="w-4 h-4" />
              </button>
            )}
          </div>
        </div>

        {/* Loading */}
        {loading && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 animate-spin text-blue-400 mx-auto mb-4" />
              <p className="text-gray-400">Loading runs...</p>
            </div>
          </div>
        )}

        {/* Error */}
        {!loading && error && (
          <div className="flex items-center justify-center py-16">
            <div className="text-center bg-gray-800 p-8 rounded-xl border border-red-700 max-w-md">
              <XCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
              <h2 className="text-lg font-bold mb-2 text-white">Error Loading Runs</h2>
              <p className="text-gray-400 text-sm mb-6">{error}</p>
              <button
                onClick={fetchRuns}
                className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
              >
                Retry
              </button>
            </div>
          </div>
        )}

        {/* Empty state */}
        {!loading && !error && filteredPermanentRuns.length === 0 && filteredRuns.length === 0 && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <FolderOpen className="w-16 h-16 text-gray-600 mx-auto mb-4" />
              <h2 className="text-xl font-bold text-white mb-2">No Runs Found</h2>
              <p className="text-gray-400 mb-6">
                {searchQuery ? 'No runs match your search.' : 'No runs have been created yet.'}
              </p>
              {searchQuery && (
                <button
                  onClick={() => setSearchQuery('')}
                  className="px-4 py-2 bg-gray-700 text-white rounded-lg hover:bg-gray-600"
                >
                  Clear Search
                </button>
              )}
            </div>
          </div>
        )}

        {/* Permanent (Main) Runs */}
        {!loading && !error && filteredPermanentRuns.length > 0 && (
          <div className="space-y-3">
            <div className="flex items-center gap-2 px-1 pt-2">
              <Database className="w-4 h-4 text-green-400" />
              <h2 className="text-sm font-semibold text-gray-200 uppercase tracking-wider">
                Permanent Runs
              </h2>
              <span className="text-xs text-gray-500">({filteredPermanentRuns.length})</span>
            </div>
            {filteredPermanentRuns.map(renderRunCard)}
          </div>
        )}

        {/* Temporary Runs */}
        {!loading && !error && filteredRuns.length > 0 && (
          <div className="space-y-3">
            <div className="flex items-center gap-2 px-1 pt-2">
              <Folder className="w-4 h-4 text-blue-400" />
              <h2 className="text-sm font-semibold text-gray-200 uppercase tracking-wider">
                Temporary Runs
              </h2>
              <span className="text-xs text-gray-500">({filteredRuns.length})</span>
            </div>
            {filteredRuns.map(renderRunCard)}
          </div>
        )}

        {/* Summary */}
        {!loading && !error && (runs.length > 0 || permanentRuns.length > 0) && (
          <div className="text-center text-sm text-gray-500">
            {permanentRuns.length} permanent, {runs.length} temporary run{runs.length !== 1 ? 's' : ''}
            {searchQuery && ` (filtered: ${filteredPermanentRuns.length} permanent, ${filteredRuns.length} temporary)`}
          </div>
        )}
      </div>
    </div>
  );
};

export default RunManagerPage;
