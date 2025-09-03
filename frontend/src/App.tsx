// src/App.tsx
import React from 'react';
import { ProteinUpload } from './components/ProteinUpload';

function App() {
  const handleMoleculesPaired = (proteinInput: string, ligandInput: string) => {
    console.log('Ready for docking:', { 
      protein: proteinInput, 
      ligand: ligandInput 
    });
    // TODO: Navigate to docking page with these inputs
  };

  return (
    <div className="min-h-screen bg-gray-50">
      <ProteinUpload onMoleculesPaired={handleMoleculesPaired} />
    </div>
  );
}

export default App;
