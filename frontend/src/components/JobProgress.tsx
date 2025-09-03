// src/components/JobProgress.tsx
import { useState, useEffect } from 'react';
import { RefreshCw, CheckCircle, XCircle, Clock } from 'lucide-react';

interface JobProgressProps {
  jobId: string;
  onComplete?: (result: any) => void;
}

interface LogEntry {
  timestamp: string;
  level: string;
  message: string;
}

export const JobProgress = ({ jobId, onComplete }: JobProgressProps) => {
  const [status, setStatus] = useState<any>(null);
  const [logs, setLogs] = useState<LogEntry[]>([]);
  const [isConnected, setIsConnected] = useState(false);

  useEffect(() => {
    // Poll job status
    const statusInterval = setInterval(async () => {
      try {
        const response = await fetch(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/jobs/${jobId}/status`);
        if (response.ok) {
          const statusData = await response.json();
          setStatus(statusData);
          
          if (statusData.status === 'completed' && onComplete) {
            onComplete(statusData);
            clearInterval(statusInterval);
          }
        }
      } catch (error) {
        console.error('Failed to fetch status:', error);
      }
    }, 1000);

    // Stream logs using Server-Sent Events
    const eventSource = new EventSource(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/jobs/${jobId}/logs`);
    
    eventSource.onopen = () => {
      setIsConnected(true);
    };
    
    eventSource.onmessage = (event) => {
      try {
        const logEntry = JSON.parse(event.data);
        setLogs(prev => [...prev, logEntry]);
      } catch (error) {
        console.error('Failed to parse log entry:', error);
      }
    };
    
    eventSource.onerror = () => {
      setIsConnected(false);
    };

    return () => {
      clearInterval(statusInterval);
      eventSource.close();
    };
  }, [jobId]);

  const getStatusIcon = () => {
    if (!status) return <Clock className="w-5 h-5 text-gray-400" />;
    
    switch (status.status) {
      case 'queued':
        return <Clock className="w-5 h-5 text-yellow-600" />;
      case 'running':
        return <RefreshCw className="w-5 h-5 text-blue-600 animate-spin" />;
      case 'completed':
        return <CheckCircle className="w-5 h-5 text-green-600" />;
      case 'failed':
        return <XCircle className="w-5 h-5 text-red-600" />;
      default:
        return <Clock className="w-5 h-5 text-gray-400" />;
    }
  };

  const getLevelColor = (level: string) => {
    switch (level) {
      case 'SUCCESS': return 'text-green-600';
      case 'ERROR': return 'text-red-600';
      case 'WARNING': return 'text-yellow-600';
      default: return 'text-gray-600';
    }
  };

  const getLevelIcon = (level: string) => {
    switch (level) {
      case 'SUCCESS': return '✅';
      case 'ERROR': return '❌';
      case 'WARNING': return '⚠️';
      default: return 'ℹ️';
    }
  };

  return (
    <div className="max-w-4xl mx-auto p-6">
      {/* Header */}
      <div className="flex items-center justify-between mb-6">
        <div className="flex items-center space-x-3">
          {getStatusIcon()}
          <div>
            <h1 className="text-2xl font-bold">Molecular Docking Progress</h1>
            <p className="text-sm text-gray-600">Job ID: {jobId}</p>
          </div>
        </div>
        <div className={`px-3 py-1 rounded-full text-sm font-medium ${
          isConnected ? 'bg-green-100 text-green-800' : 'bg-red-100 text-red-800'
        }`}>
          {isConnected ? '🔗 Connected' : '❌ Disconnected'}
        </div>
      </div>

      {/* Progress Bar */}
      {status && (
        <div className="mb-6">
          <div className="flex justify-between items-center mb-2">
            <span className="text-sm font-medium text-gray-700">
              {status.message}
            </span>
            <span className="text-sm text-gray-500">
              {status.progress}%
            </span>
          </div>
          <div className="w-full bg-gray-200 rounded-full h-3">
            <div 
              className={`h-3 rounded-full transition-all duration-300 ${
                status.status === 'completed' ? 'bg-green-600' :
                status.status === 'failed' ? 'bg-red-600' : 'bg-blue-600'
              }`}
              style={{ width: `${status.progress}%` }}
            />
          </div>
        </div>
      )}

      {/* Pipeline Logs */}
      <div className="bg-white rounded-xl border shadow-sm">
        <div className="px-6 py-4 border-b">
          <h2 className="text-lg font-semibold flex items-center space-x-2">
            <span>📋</span>
            <span>Pipeline Logs</span>
            <span className="text-sm bg-gray-100 px-2 py-1 rounded-full">
              {logs.length} entries
            </span>
          </h2>
        </div>
        <div className="p-4">
          <div className="bg-gray-900 rounded-lg p-4 h-96 overflow-y-auto font-mono text-sm">
            {logs.length === 0 ? (
              <div className="text-gray-500 text-center py-8">
                Waiting for pipeline logs...
              </div>
            ) : (
              <div className="space-y-1">
                {logs.map((log, idx) => (
                  <div key={idx} className="flex items-start space-x-2">
                    <span className="text-gray-500 shrink-0">
                      [{log.timestamp}]
                    </span>
                    <span className="shrink-0">
                      {getLevelIcon(log.level)}
                    </span>
                    <span className={`${getLevelColor(log.level)} break-all`}>
                      {log.message}
                    </span>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
      </div>

      {/* Status Info */}
      {status && (
        <div className="mt-6 grid grid-cols-1 md:grid-cols-3 gap-4">
          <div className="bg-white rounded-lg border p-4">
            <div className="text-sm text-gray-600">Status</div>
            <div className={`text-lg font-semibold capitalize ${
              status.status === 'completed' ? 'text-green-600' :
              status.status === 'failed' ? 'text-red-600' :
              status.status === 'running' ? 'text-blue-600' : 'text-yellow-600'
            }`}>
              {status.status}
            </div>
          </div>
          <div className="bg-white rounded-lg border p-4">
            <div className="text-sm text-gray-600">Progress</div>
            <div className="text-lg font-semibold">{status.progress}%</div>
          </div>
          <div className="bg-white rounded-lg border p-4">
            <div className="text-sm text-gray-600">Started</div>
            <div className="text-lg font-semibold">
              {new Date(status.created_at).toLocaleTimeString()}
            </div>
          </div>
        </div>
      )}
    </div>
  );
};
