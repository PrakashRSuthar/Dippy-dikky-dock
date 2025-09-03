// src/pages/UploadPage.tsx
import { useState } from 'react';
import { ProteinUpload } from '../components/ProteinUpload';

export const UploadPage = () => {
  const [proteinInput, setProteinInput] = useState<string>('');
  const [proteinSource, setProteinSource] = useState<'pdb' | 'upload' | 'alphafold'>('pdb');

  const handleProteinSelected = (input: string, source: 'pdb' | 'upload' | 'alphafold') => {
    setProteinInput(input);
    setProteinSource(source);
    console.log('Protein selected:', { input, source });
    // TODO: Navigate to docking page
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <ProteinUpload onProteinSelected={handleProteinSelected} />
    </div>
  );
};
