// src/pages/UploadPage.tsx
import { useState } from 'react';
import { ProteinUpload } from '../components/ProteinUpload';

export const UploadPage = () => {
  const [proteinInput, setProteinInput] = useState<string>('');
  const [ligandInput, setLigandInput] = useState<string>('');

  // Fix: ProteinUpload now expects onMoleculesPaired callback with both protein AND ligand
  const handleMoleculesPaired = (protein: string, ligand: string) => {
    setProteinInput(protein);
    setLigandInput(ligand);
    console.log('Molecules paired:', { protein, ligand });
    // TODO: Navigate to docking page
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <ProteinUpload onMoleculesPaired={handleMoleculesPaired} />
    </div>
  );
};
