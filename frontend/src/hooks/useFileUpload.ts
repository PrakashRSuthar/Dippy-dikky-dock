// src/hooks/useFileUpload.ts
import { useState } from 'react';

interface UploadResponse {
  message: string;
  file_path: string;
}

export const useFileUpload = () => {
  const [isUploading, setIsUploading] = useState(false);
  const [uploadedFileName, setUploadedFileName] = useState<string>('');

  const uploadFile = async (file: File, type: 'protein' | 'ligand'): Promise<UploadResponse | null> => {
    try {
      setIsUploading(true);
      const formData = new FormData();
      formData.append('file', file);

      const response = await fetch(`${import.meta.env.VITE_API_BASE}/api/upload/${type}`, {
        method: 'POST',
        body: formData,
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const result: UploadResponse = await response.json();
      setUploadedFileName(file.name);
      return result;
    } catch (error) {
      console.error('Upload failed:', error);
      return null;
    } finally {
      setIsUploading(false);
    }
  };

  return { uploadFile, isUploading, uploadedFileName };
};
