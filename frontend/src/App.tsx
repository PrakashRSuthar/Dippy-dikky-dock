// src/App.tsx
import React from 'react';
import { BrowserRouter, Routes, Route } from 'react-router-dom';
import Sidebar from './components/Sidebar';
import { ProteinUpload } from './components/ProteinUpload';
import { JobProgress } from './components/JobProgress';
import { ResultsPage } from './pages/ResultsPage';
import HistoryPage from './pages/HistoryPage';
import SettingsPage from './pages/SettingsPage';
import RunManagerPage from './pages/RunManagerPage';

function DockPage() {
  const [currentJobId, setCurrentJobId] = React.useState<string | null>(null);
  const [showProgress, setShowProgress] = React.useState(false);

  const handleDockingStarted = (jobId: string) => {
    setCurrentJobId(jobId);
    setShowProgress(true);
  };

  if (showProgress && currentJobId) {
    return (
      <JobProgress
        jobId={currentJobId}
        onComplete={(r: { status: string }) => {
          if (r.status === 'completed') {
            window.location.href = `/results/${currentJobId}`;
          }
        }}
      />
    );
  }

  return <ProteinUpload onDockingStarted={handleDockingStarted} />;
}

function App() {
  return (
    <BrowserRouter>
      <div className="flex min-h-screen bg-gray-950">
        <Sidebar />
        <main className="flex-1 ml-16 lg:ml-60 transition-all duration-300">
          <Routes>
            <Route path="/" element={<DockPage />} />
            <Route path="/history" element={<HistoryPage />} />
            <Route path="/runs" element={<RunManagerPage />} />
            <Route path="/settings" element={<SettingsPage />} />
            <Route path="/results/:jobId" element={<ResultsPageWrapper />} />
            <Route path="/progress/:jobId" element={<ProgressWrapper />} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  );
}

function ResultsPageWrapper() {
  const jobId = window.location.pathname.split('/results/')[1];
  return <ResultsPage jobId={jobId} onBack={() => window.location.href = '/history'} />;
}

function ProgressWrapper() {
  const jobId = window.location.pathname.split('/progress/')[1];
  return (
    <div className="min-h-screen bg-gray-950 p-6">
      <JobProgress
        jobId={jobId}
        onComplete={(r) => {
          if (r.status === 'completed') {
            window.location.href = `/results/${jobId}`;
          }
        }}
      />
    </div>
  );
}

export default App;
