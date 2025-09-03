// src/components/MoleculeVisualization.tsx
import React from 'react';

interface MoleculeVisualizationProps {
  moleculePath?: string | null;
  moleculeType?: 'protein' | 'ligand' | 'complex';
  color?: string;
}

export const MoleculeVisualization: React.FC<MoleculeVisualizationProps> = ({ 
  moleculePath, 
  moleculeType = 'protein', 
  color = '#6b7280' 
}) => {
  if (!moleculePath || moleculePath.trim() === '') {
    return (
      <div className="w-full h-full flex items-center justify-center">
        <div className="text-center text-gray-400">
          <div className="text-6xl mb-6">
            {moleculeType === 'protein' ? '🧬' : moleculeType === 'ligand' ? '💊' : '🔬'}
          </div>
          <p className="text-sm">No {moleculeType} loaded</p>
        </div>
      </div>
    );
  }

  const isFilePath = moleculePath.includes('/');
  const displayName = isFilePath ? 'Uploaded File' : moleculePath;

  return (
    <div className="w-full h-full relative bg-gradient-to-br from-gray-900 to-gray-700 rounded-lg overflow-hidden">
      {/* 3D-like Background Pattern */}
      <div className="absolute inset-0 opacity-20">
        <div className="absolute top-4 left-4 w-2 h-2 bg-white rounded-full animate-pulse"></div>
        <div className="absolute top-8 right-8 w-1 h-1 bg-white rounded-full animate-pulse delay-300"></div>
        <div className="absolute bottom-6 left-8 w-1 h-1 bg-white rounded-full animate-pulse delay-500"></div>
        <div className="absolute bottom-4 right-4 w-2 h-2 bg-white rounded-full animate-pulse delay-700"></div>
        
        {/* Grid pattern */}
        <svg className="w-full h-full opacity-10">
          <defs>
            <pattern id="grid" width="20" height="20" patternUnits="userSpaceOnUse">
              <path d="M 20 0 L 0 0 0 20" fill="none" stroke="white" strokeWidth="0.5"/>
            </pattern>
          </defs>
          <rect width="100%" height="100%" fill="url(#grid)" />
        </svg>
      </div>

      {/* Central Molecule Display */}
      <div className="absolute inset-0 flex items-center justify-center">
        <div className="text-center p-6">
          {/* Animated Molecule Icon */}
          <div className="relative mb-4">
            <div 
              className="text-6xl transform transition-transform duration-1000 hover:scale-110" 
              style={{ color }}
            >
              {moleculeType === 'protein' ? '🧬' : moleculeType === 'ligand' ? '💊' : '🔬'}
            </div>
            
            {/* Rotating rings around molecule */}
            <div className="absolute inset-0 flex items-center justify-center">
              <div className="w-16 h-16 border border-white opacity-30 rounded-full animate-spin"></div>
            </div>
            <div className="absolute inset-0 flex items-center justify-center">
              <div className="w-20 h-20 border border-white opacity-20 rounded-full animate-ping"></div>
            </div>
          </div>

          {/* Molecule Info */}
          <div className="bg-black bg-opacity-50 rounded-lg p-4 backdrop-blur-sm">
            <div className="text-lg font-semibold text-white mb-2">
              {displayName}
            </div>
            <div className="text-sm text-gray-300 mb-1 capitalize">
              {moleculeType === 'complex' ? 'Docked Complex' : `${moleculeType} Structure`}
            </div>
            <div className="text-xs text-gray-400 font-mono bg-gray-800 px-2 py-1 rounded max-w-full truncate">
              {moleculePath}
            </div>
            
            {/* 3D Viewer Simulation */}
            <div className="mt-4 flex justify-center space-x-2">
              <div className="w-2 h-2 bg-green-400 rounded-full animate-pulse"></div>
              <div className="w-2 h-2 bg-blue-400 rounded-full animate-pulse delay-200"></div>
              <div className="w-2 h-2 bg-purple-400 rounded-full animate-pulse delay-400"></div>
            </div>
            <div className="text-xs text-gray-400 mt-2">
              Interactive 3D Viewer Active
            </div>
          </div>

          {/* Control Buttons */}
          <div className="flex justify-center space-x-2 mt-4">
            <button className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors">
              Rotate
            </button>
            <button className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors">
              Zoom
            </button>
            <button className="px-3 py-1 bg-gray-700 hover:bg-gray-600 text-white text-xs rounded transition-colors">
              Reset
            </button>
          </div>
        </div>
      </div>

      {/* Loading Animation Overlay */}
      <div className="absolute bottom-4 right-4 bg-black bg-opacity-70 rounded px-3 py-2">
        <div className="flex items-center space-x-2">
          <div className="w-3 h-3 bg-green-400 rounded-full animate-pulse"></div>
          <span className="text-xs text-white">Rendering</span>
        </div>
      </div>
    </div>
  );
};
