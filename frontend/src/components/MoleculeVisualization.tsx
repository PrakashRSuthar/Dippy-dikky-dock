// src/components/MoleculeVisualization.tsx
import React, { Suspense } from 'react';
import { Canvas } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';

interface MoleculeVisualizationProps {
  proteinPath?: string | null;
  ligandPath?: string | null;
}

// Simple placeholder for molecular visualization
function MoleculeViewer({ modelPath, color = '#4f46e5' }: { modelPath: string; color?: string }) {
  return (
    <mesh>
      <sphereGeometry args={[1, 32, 32]} />
      <meshStandardMaterial color={color} />
    </mesh>
  );
}

export function MoleculeVisualization({ proteinPath, ligandPath }: MoleculeVisualizationProps) {
  return (
    <Canvas camera={{ position: [0, 0, 5], fov: 45 }} className="w-full h-full">
      <ambientLight intensity={0.5} />
      <directionalLight position={[5, 5, 5]} intensity={1} />
      <Environment preset="studio" />
      
      <Suspense fallback={null}>
        {proteinPath && (
          <group position={[-1.5, 0, 0]}>
            <MoleculeViewer modelPath={proteinPath} color="#3b82f6" />
          </group>
        )}
        {ligandPath && (
          <group position={[1.5, 0, 0]}>
            <MoleculeViewer modelPath={ligandPath} color="#10b981" />
          </group>
        )}
      </Suspense>
      
      <OrbitControls enableZoom enablePan enableRotate />
    </Canvas>
  );
}
