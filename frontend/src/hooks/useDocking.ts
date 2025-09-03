// src/hooks/useDocking.ts
import { useState } from 'react';

interface DockingRequest {
  protein_input: string;
  ligand_input: string;
  job_name?: string;
  cleaning_policy?: {
    keep_waters: boolean;
    keep_ions: boolean;
    keep_solvents: boolean;
    keep_cofactors: boolean;
  };
}

interface DockingResponse {
  job_id: string;
  message: string;
  status: string;
}

export const useDocking = () => {
  const [isStarting, setIsStarting] = useState(false);

  const startDocking = async (request: DockingRequest): Promise<string | null> => {
    try {
      setIsStarting(true);
      const response = await fetch(`${import.meta.env.VITE_API_BASE || 'http://localhost:8000'}/api/dock`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(request),
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const result: DockingResponse = await response.json();
      console.log('Docking started:', result);
      return result.job_id;
    } catch (error) {
      console.error('Failed to start docking:', error);
      alert('Failed to start docking. Please check the backend connection.');
      return null;
    } finally {
      setIsStarting(false);
    }
  };

  return { startDocking, isStarting };
};
