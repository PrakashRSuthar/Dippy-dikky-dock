// src/hooks/useJobStatus.ts
import { useState, useEffect } from 'react';

interface JobStatus {
  job_id: string;
  status: 'queued' | 'running' | 'completed' | 'failed';
  progress: number;
  message: string;
  created_at: string;
  completed_at?: string;
}

export const useJobStatus = (jobId: string) => {
  const [jobStatus, setJobStatus] = useState<JobStatus | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const fetchStatus = async () => {
      try {
        const response = await fetch(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/jobs/${jobId}/status`);
        if (response.ok) {
          const status = await response.json();
          setJobStatus(status);
        }
      } catch (error) {
        console.error('Failed to fetch job status:', error);
      } finally {
        setLoading(false);
      }
    };

    fetchStatus();
    
    // Poll for status updates if job is not completed
    const interval = setInterval(() => {
      if (jobStatus?.status === 'running' || jobStatus?.status === 'queued') {
        fetchStatus();
      }
    }, 2000);

    return () => clearInterval(interval);
  }, [jobId, jobStatus?.status]);

  return { jobStatus, loading };
};
