// src/hooks/useApiHealth.ts
import { useState, useEffect } from 'react';

interface HealthResponse {
  message: string;
  version: string;
  status: string;
}

export const useApiHealth = () => {
  const [isHealthy, setIsHealthy] = useState(false);
  const [isChecking, setIsChecking] = useState(true);

  const checkHealth = async () => {
    try {
      setIsChecking(true);
      const response = await fetch(`${import.meta.env.VITE_API_BASE}/`);
      
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      
      const data: HealthResponse = await response.json();
      setIsHealthy(data.status === 'running');
    } catch (error) {
      console.error('Health check failed:', error);
      setIsHealthy(false);
    } finally {
      setIsChecking(false);
    }
  };

  useEffect(() => {
    checkHealth();
    const interval = setInterval(checkHealth, 30000); // Check every 30 seconds
    return () => clearInterval(interval);
  }, []);

  return { isHealthy, isChecking };
};
