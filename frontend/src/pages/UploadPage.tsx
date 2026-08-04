// UploadPage.tsx — Fixed: prop now matches ProteinUpload.onDockingStarted(jobId: string)
import { useState } from 'react';
import { ProteinUpload } from '../components/ProteinUpload';

export const UploadPage = () => {
  const [_activeJobId, _setActiveJobId] = useState<string | null>(null);

  const handleDockingStarted = (jobId: string) => {
    _setActiveJobId(jobId);
    console.log('Job started:', jobId);
    // Navigate to progress / results page here
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <ProteinUpload onDockingStarted={handleDockingStarted} />
    </div>
  );
};