// src/App.tsx
import React, { useState } from 'react';
import { ProteinUpload } from './components/ProteinUpload';
import { JobProgress } from './components/JobProgress';
import { ResultsPage } from './pages/ResultsPage';

type AppState = 'upload' | 'progress' | 'results';

function App() {
  const [currentState, setCurrentState] = useState<AppState>('upload');
  const [currentJobId, setCurrentJobId] = useState<string | null>(null);

  const handleDockingStarted = (jobId: string) => {
    setCurrentJobId(jobId);
    setCurrentState('progress');
  };

  const handleJobComplete = (result: any) => {
    if (result.status === 'completed') {
      setCurrentState('results');
    }
  };

  const handleBackToUpload = () => {
    setCurrentJobId(null);
    setCurrentState('upload');
  };

  if (currentState === 'progress' && currentJobId) {
    return (
      <div className="min-h-screen bg-gray-50">
        <JobProgress 
          jobId={currentJobId} 
          onComplete={handleJobComplete}
        />
      </div>
    );
  }

  if (currentState === 'results' && currentJobId) {
    return <ResultsPage jobId={currentJobId} onBack={handleBackToUpload} />;
  }

  return (
    <div className="min-h-screen bg-gray-50">
      <ProteinUpload onDockingStarted={handleDockingStarted} />
    </div>
  );
}

export default App;
