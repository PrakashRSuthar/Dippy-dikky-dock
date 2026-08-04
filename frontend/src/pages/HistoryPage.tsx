import { useState, useEffect, useCallback, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import {
  Search, RefreshCw, Trash2, Edit3, Eye, ChevronRight,
  CheckCircle2, XCircle, Circle, Clock, BarChart3, TrendingUp, Activity,
  FolderOpen, Filter, X
} from 'lucide-react';

const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

type JobStatus = 'queued' | 'running' | 'completed' | 'failed';

interface HistoryJob {
  id: string;
  name: string;
  status: JobStatus;
  protein_name?: string;
  ligand_name?: string;
  best_affinity?: number;
  rmsd?: number;
  duration?: string;
  created_at: string;
  updated_at?: string;
}

interface Stats {
  total_jobs: number;
  success_rate: number;
  avg_rmsd: number;
  completed_count: number;
  failed_count: number;
  running_count: number;
}

const statusConfig: Record<JobStatus, { color: string; bg: string; icon: React.ReactNode }> = {
  completed: {
    color: 'text-green-400',
    bg: 'bg-green-900/30 border-green-700',
    icon: <CheckCircle2 className="w-3.5 h-3.5" />
  },
  failed: {
    color: 'text-red-400',
    bg: 'bg-red-900/30 border-red-700',
    icon: <XCircle className="w-3.5 h-3.5" />
  },
  running: {
    color: 'text-yellow-400',
    bg: 'bg-yellow-900/30 border-yellow-700',
    icon: <Circle className="w-3.5 h-3.5 animate-pulse" />
  },
  queued: {
    color: 'text-gray-400',
    bg: 'bg-gray-800 border-gray-600',
    icon: <Clock className="w-3.5 h-3.5" />
  }
};

function StatusBadge({ status }: { status: JobStatus }) {
  const config = statusConfig[status];
  return (
    <span className={`inline-flex items-center gap-1.5 px-2.5 py-1 rounded-full text-xs font-medium border ${config.bg} ${config.color}`}>
      {config.icon}
      {status.charAt(0).toUpperCase() + status.slice(1)}
    </span>
  );
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

export const HistoryPage = () => {
  const navigate = useNavigate();

  const [jobs, setJobs] = useState<HistoryJob[]>([]);
  const [stats, setStats] = useState<Stats | null>(null);
  const [loading, setLoading] = useState(true);
  const [loadingMore, setLoadingMore] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [offset, setOffset] = useState(0);
  const [hasMore, setHasMore] = useState(true);
  const [searchQuery, setSearchQuery] = useState('');
  const [statusFilter, setStatusFilter] = useState<JobStatus | 'all'>('all');
  const [editingId, setEditingId] = useState<string | null>(null);
  const [editingName, setEditingName] = useState('');
  const [deleteConfirmId, setDeleteConfirmId] = useState<string | null>(null);
  const [toast, setToast] = useState<{ message: string; type: 'success' | 'error' } | null>(null);

  const limit = 50;

  const showToast = useCallback((message: string, type: 'success' | 'error' = 'success') => {
    setToast({ message, type });
    setTimeout(() => setToast(null), 3000);
  }, []);

  const fetchStats = useCallback(async () => {
    try {
      const res = await fetch(`${apiBase}/api/stats`);
      if (res.ok) {
        const data = await res.json();
        setStats(data);
      }
    } catch (e) {
      console.error('Failed to fetch stats:', e);
    }
  }, []);

  const fetchJobs = useCallback(async (currentOffset: number, append: boolean = false) => {
    try {
      const res = await fetch(`${apiBase}/api/history?limit=${limit}&offset=${currentOffset}`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      const newJobs = data.jobs || data || [];
      
      if (append) {
        setJobs(prev => [...prev, ...newJobs]);
      } else {
        setJobs(newJobs);
      }
      
      setHasMore(newJobs.length === limit);
    } catch (e: any) {
      setError(e?.message || 'Failed to load history');
    }
  }, []);

  const loadInitial = useCallback(async () => {
    setLoading(true);
    setError(null);
    setOffset(0);
    await Promise.all([fetchJobs(0), fetchStats()]);
    setLoading(false);
  }, [fetchJobs, fetchStats]);

  useEffect(() => {
    loadInitial();
  }, [loadInitial]);

  const handleLoadMore = useCallback(async () => {
    const newOffset = offset + limit;
    setLoadingMore(true);
    await fetchJobs(newOffset, true);
    setOffset(newOffset);
    setLoadingMore(false);
  }, [offset, fetchJobs]);

  const handleRefresh = useCallback(() => {
    loadInitial();
  }, [loadInitial]);

  const handleView = useCallback((jobId: string) => {
    navigate(`/results/${jobId}`);
  }, [navigate]);

  const handleStartEdit = useCallback((job: HistoryJob) => {
    setEditingId(job.id);
    setEditingName(job.name);
  }, []);

  const handleCancelEdit = useCallback(() => {
    setEditingId(null);
    setEditingName('');
  }, []);

  const handleSaveName = useCallback(async (jobId: string) => {
    try {
      const res = await fetch(`${apiBase}/api/jobs/${jobId}`, {
        method: 'PATCH',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: editingName })
      });
      if (!res.ok) throw new Error('Failed to rename');
      setJobs(prev => prev.map(j => j.id === jobId ? { ...j, name: editingName } : j));
      setEditingId(null);
      setEditingName('');
      showToast('Job renamed successfully');
    } catch (e) {
      showToast('Failed to rename job', 'error');
    }
  }, [editingName, showToast]);

  const handleDelete = useCallback(async (jobId: string) => {
    try {
      const res = await fetch(`${apiBase}/api/jobs/${jobId}`, { method: 'DELETE' });
      if (!res.ok) throw new Error('Failed to delete');
      setJobs(prev => prev.filter(j => j.id !== jobId));
      setDeleteConfirmId(null);
      showToast('Job deleted successfully');
      fetchStats();
    } catch (e) {
      showToast('Failed to delete job', 'error');
    }
  }, [showToast, fetchStats]);

  const filteredJobs = useMemo(() => {
    return jobs.filter(job => {
      const matchesSearch = searchQuery === '' || 
        job.name.toLowerCase().includes(searchQuery.toLowerCase()) ||
        job.id.toLowerCase().includes(searchQuery.toLowerCase()) ||
        (job.protein_name?.toLowerCase().includes(searchQuery.toLowerCase()) ?? false) ||
        (job.ligand_name?.toLowerCase().includes(searchQuery.toLowerCase()) ?? false);
      const matchesStatus = statusFilter === 'all' || job.status === statusFilter;
      return matchesSearch && matchesStatus;
    });
  }, [jobs, searchQuery, statusFilter]);

  return (
    <div className="min-h-screen bg-gray-900 text-gray-100">
      {/* Toast */}
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
              <Activity className="w-6 h-6 text-blue-400" />
              <h1 className="text-2xl font-bold">Job History</h1>
            </div>
            <button
              onClick={handleRefresh}
              disabled={loading}
              className="flex items-center gap-2 px-4 py-2 bg-gray-700 hover:bg-gray-600 rounded-lg transition-colors disabled:opacity-50"
            >
              <RefreshCw className={`w-4 h-4 ${loading ? 'animate-spin' : ''}`} />
              Refresh
            </button>
          </div>
        </div>
      </div>

      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-6">
        {/* Stats Summary */}
        {stats && (
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            <div className="bg-gray-800 rounded-xl border border-gray-700 p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-gray-400 mb-1">Total Jobs</p>
                  <p className="text-2xl font-bold text-white">{stats.total_jobs}</p>
                </div>
                <BarChart3 className="w-8 h-8 text-blue-500 opacity-60" />
              </div>
            </div>
            <div className="bg-gray-800 rounded-xl border border-gray-700 p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-gray-400 mb-1">Success Rate</p>
                  <p className="text-2xl font-bold text-green-400">{stats.success_rate.toFixed(1)}%</p>
                </div>
                <TrendingUp className="w-8 h-8 text-green-500 opacity-60" />
              </div>
            </div>
            <div className="bg-gray-800 rounded-xl border border-gray-700 p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-gray-400 mb-1">Avg RMSD</p>
                  <p className="text-2xl font-bold text-purple-400">{stats.avg_rmsd.toFixed(2)} Å</p>
                </div>
                <Activity className="w-8 h-8 text-purple-500 opacity-60" />
              </div>
            </div>
            <div className="bg-gray-800 rounded-xl border border-gray-700 p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-gray-400 mb-1">Completed</p>
                  <p className="text-2xl font-bold text-blue-400">{stats.completed_count}</p>
                </div>
                <CheckCircle2 className="w-8 h-8 text-blue-500 opacity-60" />
              </div>
            </div>
          </div>
        )}

        {/* Search and Filter */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-4 mb-6">
          <div className="flex flex-col sm:flex-row gap-4">
            <div className="flex-1 relative">
              <Search className="absolute left-3 top-1/2 transform -translate-y-1/2 w-4 h-4 text-gray-400" />
              <input
                type="text"
                placeholder="Search by name, ID, protein, or ligand..."
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
            <div className="flex items-center gap-2">
              <Filter className="w-4 h-4 text-gray-400" />
              <select
                value={statusFilter}
                onChange={(e) => setStatusFilter(e.target.value as JobStatus | 'all')}
                className="px-4 py-2 bg-gray-700 border border-gray-600 rounded-lg text-white focus:outline-none focus:ring-2 focus:ring-blue-500"
              >
                <option value="all">All Status</option>
                <option value="completed">Completed</option>
                <option value="failed">Failed</option>
                <option value="running">Running</option>
                <option value="queued">Queued</option>
              </select>
            </div>
          </div>
        </div>

        {/* Loading State */}
        {loading && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <RefreshCw className="w-8 h-8 animate-spin text-blue-400 mx-auto mb-4" />
              <p className="text-gray-400">Loading job history...</p>
            </div>
          </div>
        )}

        {/* Error State */}
        {!loading && error && (
          <div className="flex items-center justify-center py-16">
            <div className="text-center bg-gray-800 p-8 rounded-xl border border-red-700 max-w-md">
              <XCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
              <h2 className="text-lg font-bold mb-2 text-white">Error Loading History</h2>
              <p className="text-gray-400 text-sm mb-6">{error}</p>
              <button
                onClick={handleRefresh}
                className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
              >
                Retry
              </button>
            </div>
          </div>
        )}

        {/* Empty State */}
        {!loading && !error && filteredJobs.length === 0 && (
          <div className="flex items-center justify-center py-24">
            <div className="text-center">
              <FolderOpen className="w-16 h-16 text-gray-600 mx-auto mb-4" />
              <h2 className="text-xl font-bold text-white mb-2">No Jobs Found</h2>
              <p className="text-gray-400 mb-6">
                {searchQuery || statusFilter !== 'all'
                  ? 'No jobs match your search criteria.'
                  : 'Start your first docking job to see results here.'}
              </p>
              {(searchQuery || statusFilter !== 'all') && (
                <button
                  onClick={() => { setSearchQuery(''); setStatusFilter('all'); }}
                  className="px-4 py-2 bg-gray-700 text-white rounded-lg hover:bg-gray-600"
                >
                  Clear Filters
                </button>
              )}
            </div>
          </div>
        )}

        {/* Jobs Table */}
        {!loading && !error && filteredJobs.length > 0 && (
          <>
            <div className="bg-gray-800 rounded-xl border border-gray-700 overflow-hidden">
              <div className="overflow-x-auto">
                <table className="w-full">
                  <thead>
                    <tr className="bg-gray-750 border-b border-gray-700">
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Name</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Status</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Protein</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Ligand</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Best Affinity</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">RMSD</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Duration</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Date</th>
                      <th className="px-4 py-3 text-left text-xs font-semibold text-gray-300 uppercase tracking-wider">Actions</th>
                    </tr>
                  </thead>
                  <tbody className="divide-y divide-gray-700">
                    {filteredJobs.map((job) => (
                      <tr key={job.id} className="hover:bg-gray-750 transition-colors">
                        <td className="px-4 py-3">
                          {editingId === job.id ? (
                            <div className="flex items-center gap-2">
                              <input
                                type="text"
                                value={editingName}
                                onChange={(e) => setEditingName(e.target.value)}
                                onKeyDown={(e) => {
                                  if (e.key === 'Enter') handleSaveName(job.id);
                                  if (e.key === 'Escape') handleCancelEdit();
                                }}
                                className="flex-1 px-2 py-1 bg-gray-700 border border-gray-600 rounded text-white text-sm focus:outline-none focus:ring-2 focus:ring-blue-500"
                                autoFocus
                              />
                              <button onClick={() => handleSaveName(job.id)} className="text-green-400 hover:text-green-300">
                                <CheckCircle2 className="w-4 h-4" />
                              </button>
                              <button onClick={handleCancelEdit} className="text-gray-400 hover:text-white">
                                <X className="w-4 h-4" />
                              </button>
                            </div>
                          ) : (
                            <span className="font-medium text-white">{job.name || 'Untitled Job'}</span>
                          )}
                        </td>
                        <td className="px-4 py-3">
                          <StatusBadge status={job.status} />
                        </td>
                        <td className="px-4 py-3 text-sm text-gray-300">{job.protein_name || '—'}</td>
                        <td className="px-4 py-3 text-sm text-gray-300">{job.ligand_name || '—'}</td>
                        <td className="px-4 py-3 text-sm font-mono text-green-400">
                          {job.best_affinity != null ? `${job.best_affinity.toFixed(2)} kcal/mol` : '—'}
                        </td>
                        <td className="px-4 py-3 text-sm font-mono text-purple-400">
                          {job.rmsd != null ? `${job.rmsd.toFixed(2)} Å` : '—'}
                        </td>
                        <td className="px-4 py-3 text-sm text-gray-300">{job.duration || '—'}</td>
                        <td className="px-4 py-3 text-sm text-gray-400">{formatDate(job.created_at)}</td>
                        <td className="px-4 py-3">
                          {deleteConfirmId === job.id ? (
                            <div className="flex items-center gap-2">
                              <span className="text-xs text-red-400">Delete?</span>
                              <button
                                onClick={() => handleDelete(job.id)}
                                className="px-2 py-1 bg-red-600 text-white text-xs rounded hover:bg-red-700"
                              >
                                Yes
                              </button>
                              <button
                                onClick={() => setDeleteConfirmId(null)}
                                className="px-2 py-1 bg-gray-600 text-white text-xs rounded hover:bg-gray-500"
                              >
                                No
                              </button>
                            </div>
                          ) : (
                            <div className="flex items-center gap-1">
                              <button
                                onClick={() => handleView(job.id)}
                                className="p-1.5 text-gray-400 hover:text-blue-400 hover:bg-gray-700 rounded transition-colors"
                                title="View Results"
                              >
                                <Eye className="w-4 h-4" />
                              </button>
                              <button
                                onClick={() => handleStartEdit(job)}
                                className="p-1.5 text-gray-400 hover:text-yellow-400 hover:bg-gray-700 rounded transition-colors"
                                title="Rename"
                              >
                                <Edit3 className="w-4 h-4" />
                              </button>
                              <button
                                onClick={() => setDeleteConfirmId(job.id)}
                                className="p-1.5 text-gray-400 hover:text-red-400 hover:bg-gray-700 rounded transition-colors"
                                title="Delete"
                              >
                                <Trash2 className="w-4 h-4" />
                              </button>
                            </div>
                          )}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>

            {/* Load More */}
            {hasMore && (
              <div className="flex justify-center mt-6">
                <button
                  onClick={handleLoadMore}
                  disabled={loadingMore}
                  className="flex items-center gap-2 px-6 py-3 bg-gray-800 border border-gray-700 text-white rounded-lg hover:bg-gray-700 transition-colors disabled:opacity-50"
                >
                  {loadingMore ? (
                    <RefreshCw className="w-4 h-4 animate-spin" />
                  ) : (
                    <ChevronRight className="w-4 h-4" />
                  )}
                  {loadingMore ? 'Loading...' : 'Load More'}
                </button>
              </div>
            )}

            {/* Results count */}
            <div className="text-center mt-4 text-sm text-gray-500">
              Showing {filteredJobs.length} of {jobs.length} jobs
            </div>
          </>
        )}
      </div>
    </div>
  );
};

export default HistoryPage;
