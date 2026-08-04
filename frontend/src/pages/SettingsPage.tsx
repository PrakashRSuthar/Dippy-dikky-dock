import { useState, useEffect, useCallback } from 'react';
import {
  Settings, Save, RefreshCw, Box, Droplets, Clock, Palette, Info,
  Zap, XCircle, ExternalLink
} from 'lucide-react';

const apiBase = (import.meta as any).env?.VITE_API_BASE || 'http://localhost:8000';

interface DockingSettings {
  exhaustiveness: number;
  num_modes: number;
  energy_range: number;
}

interface BoxSettings {
  auto_box: boolean;
  min_axis: number;
  max_axis: number;
}

interface CleaningSettings {
  keep_waters: boolean;
  keep_ions: boolean;
  keep_solvents: boolean;
  keep_cofactors: boolean;
}

interface RetentionSettings {
  retention: 'permanent' | '7days' | 'immediate';
}

interface AppSettings {
  docking: DockingSettings;
  box: BoxSettings;
  cleaning: CleaningSettings;
  retention: RetentionSettings;
}

const defaultSettings: AppSettings = {
  docking: {
    exhaustiveness: 8,
    num_modes: 9,
    energy_range: 3
  },
  box: {
    auto_box: true,
    min_axis: 20,
    max_axis: 30
  },
  cleaning: {
    keep_waters: false,
    keep_ions: false,
    keep_solvents: false,
    keep_cofactors: true
  },
  retention: {
    retention: 'permanent'
  }
};

interface SliderProps {
  label: string;
  value: number;
  min: number;
  max: number;
  onChange: (value: number) => void;
  description?: string;
}

function Slider({ label, value, min, max, onChange, description }: SliderProps) {
  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <label className="text-sm font-medium text-gray-200">{label}</label>
        <span className="text-sm font-mono text-blue-400 bg-gray-700 px-2 py-0.5 rounded">{value}</span>
      </div>
      {description && <p className="text-xs text-gray-400">{description}</p>}
      <input
        type="range"
        min={min}
        max={max}
        value={value}
        onChange={(e) => onChange(parseInt(e.target.value))}
        className="w-full h-2 bg-gray-700 rounded-lg appearance-none cursor-pointer accent-blue-500"
      />
      <div className="flex justify-between text-xs text-gray-500">
        <span>{min}</span>
        <span>{max}</span>
      </div>
    </div>
  );
}

interface ToggleProps {
  label: string;
  checked: boolean;
  onChange: (checked: boolean) => void;
  description?: string;
}

function Toggle({ label, checked, onChange, description }: ToggleProps) {
  return (
    <div className="flex items-center justify-between py-2">
      <div className="flex-1">
        <label className="text-sm font-medium text-gray-200">{label}</label>
        {description && <p className="text-xs text-gray-400 mt-0.5">{description}</p>}
      </div>
      <button
        type="button"
        onClick={() => onChange(!checked)}
        className={`relative inline-flex h-6 w-11 items-center rounded-full transition-colors ${
          checked ? 'bg-blue-600' : 'bg-gray-600'
        }`}
      >
        <span
          className={`inline-block h-4 w-4 transform rounded-full bg-white transition-transform ${
            checked ? 'translate-x-6' : 'translate-x-1'
          }`}
        />
      </button>
    </div>
  );
}

interface RadioGroupProps {
  label: string;
  options: { value: string; label: string; description?: string }[];
  selected: string;
  onChange: (value: string) => void;
}

function RadioGroup({ label, options, selected, onChange }: RadioGroupProps) {
  return (
    <div className="space-y-3">
      <label className="text-sm font-medium text-gray-200">{label}</label>
      <div className="space-y-2">
        {options.map((option) => (
          <label
            key={option.value}
            className={`flex items-start gap-3 p-3 rounded-lg border cursor-pointer transition-colors ${
              selected === option.value
                ? 'bg-blue-900/30 border-blue-600'
                : 'bg-gray-800 border-gray-700 hover:border-gray-600'
            }`}
          >
            <input
              type="radio"
              name={label}
              value={option.value}
              checked={selected === option.value}
              onChange={() => onChange(option.value)}
              className="mt-0.5 accent-blue-500"
            />
            <div>
              <span className="text-sm font-medium text-white">{option.label}</span>
              {option.description && (
                <p className="text-xs text-gray-400 mt-0.5">{option.description}</p>
              )}
            </div>
          </label>
        ))}
      </div>
    </div>
  );
}

export const SettingsPage = () => {
  const [settings, setSettings] = useState<AppSettings>(defaultSettings);
  const [loading, setLoading] = useState(true);
  const [saving, setSaving] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [toast, setToast] = useState<{ message: string; type: 'success' | 'error' } | null>(null);

  const showToast = useCallback((message: string, type: 'success' | 'error' = 'success') => {
    setToast({ message, type });
    setTimeout(() => setToast(null), 3000);
  }, []);

  const fetchSettings = useCallback(async () => {
    try {
      setLoading(true);
      setError(null);
      const res = await fetch(`${apiBase}/api/settings`);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const data = await res.json();
      const s = data.settings || data;
      
      // Merge with defaults to ensure all fields exist
      setSettings({
        docking: { ...defaultSettings.docking, ...s.docking },
        box: { ...defaultSettings.box, ...s.box },
        cleaning: { ...defaultSettings.cleaning, ...s.cleaning },
        retention: { ...defaultSettings.retention, ...s.retention }
      });
    } catch (e: any) {
      setError(e?.message || 'Failed to load settings');
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    fetchSettings();
  }, [fetchSettings]);

  const handleSave = useCallback(async () => {
    try {
      setSaving(true);
      const res = await fetch(`${apiBase}/api/settings`, {
        method: 'PUT',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ settings })
      });
      if (!res.ok) throw new Error('Failed to save settings');
      showToast('Settings saved successfully');
    } catch (e: any) {
      showToast(e?.message || 'Failed to save settings', 'error');
    } finally {
      setSaving(false);
    }
  }, [settings, showToast]);

  const updateDocking = useCallback((key: keyof DockingSettings, value: number) => {
    setSettings(prev => ({
      ...prev,
      docking: { ...prev.docking, [key]: value }
    }));
  }, []);

  const updateBox = useCallback((key: keyof BoxSettings, value: boolean | number) => {
    setSettings(prev => ({
      ...prev,
      box: { ...prev.box, [key]: value }
    }));
  }, []);

  const updateCleaning = useCallback((key: keyof CleaningSettings, value: boolean) => {
    setSettings(prev => ({
      ...prev,
      cleaning: { ...prev.cleaning, [key]: value }
    }));
  }, []);

  const updateRetention = useCallback((value: string) => {
    setSettings(prev => ({
      ...prev,
      retention: { retention: value as AppSettings['retention']['retention'] }
    }));
  }, []);

  if (loading) {
    return (
      <div className="min-h-screen bg-gray-900 flex items-center justify-center">
        <div className="text-center">
          <RefreshCw className="w-8 h-8 animate-spin text-blue-400 mx-auto mb-4" />
          <p className="text-gray-400">Loading settings...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="min-h-screen bg-gray-900 flex items-center justify-center">
        <div className="text-center bg-gray-800 p-8 rounded-xl border border-red-700 max-w-md">
          <XCircle className="w-12 h-12 text-red-500 mx-auto mb-4" />
          <h2 className="text-lg font-bold mb-2 text-white">Error Loading Settings</h2>
          <p className="text-gray-400 text-sm mb-6">{error}</p>
          <button
            onClick={fetchSettings}
            className="px-4 py-2 bg-blue-600 text-white rounded-lg hover:bg-blue-700"
          >
            Retry
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-gray-900 text-gray-100">
      {/* Toast */}
      {toast && (
        <div className={`fixed top-4 right-4 z-50 px-4 py-3 rounded-lg shadow-lg transition-all ${
          toast.type === 'success' ? 'bg-green-600 text-white' : 'bg-red-600 text-white'
        }`}>
          {toast.message}
        </div>
      )}

      {/* Header */}
      <div className="bg-gray-800 border-b border-gray-700 sticky top-0 z-10">
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center gap-3">
              <Settings className="w-6 h-6 text-blue-400" />
              <h1 className="text-2xl font-bold">Settings</h1>
            </div>
            <button
              onClick={handleSave}
              disabled={saving}
              className="flex items-center gap-2 px-4 py-2 bg-blue-600 hover:bg-blue-700 text-white rounded-lg transition-colors disabled:opacity-50"
            >
              {saving ? (
                <RefreshCw className="w-4 h-4 animate-spin" />
              ) : (
                <Save className="w-4 h-4" />
              )}
              {saving ? 'Saving...' : 'Save Changes'}
            </button>
          </div>
        </div>
      </div>

      <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-8 space-y-6">
        {/* Docking Engine */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-blue-900/50 rounded-lg">
              <Zap className="w-5 h-5 text-blue-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">Docking Engine</h2>
              <p className="text-sm text-gray-400">Configure AutoDock Vina parameters</p>
            </div>
          </div>
          
          <div className="space-y-6">
            <Slider
              label="Exhaustiveness"
              value={settings.docking.exhaustiveness}
              min={1}
              max={32}
              onChange={(v) => updateDocking('exhaustiveness', v)}
              description="Higher values search more thoroughly but take longer. Default: 8"
            />
            <Slider
              label="Num Modes"
              value={settings.docking.num_modes}
              min={1}
              max={20}
              onChange={(v) => updateDocking('num_modes', v)}
              description="Maximum number of binding modes to generate. Default: 9"
            />
            <Slider
              label="Energy Range"
              value={settings.docking.energy_range}
              min={1}
              max={10}
              onChange={(v) => updateDocking('energy_range', v)}
              description="Maximum energy difference from best mode (kcal/mol). Default: 3"
            />
          </div>
        </div>

        {/* Box Settings */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-purple-900/50 rounded-lg">
              <Box className="w-5 h-5 text-purple-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">Box Settings</h2>
              <p className="text-sm text-gray-400">Configure the docking search space</p>
            </div>
          </div>
          
          <div className="space-y-6">
            <Toggle
              label="Auto Box"
              checked={settings.box.auto_box}
              onChange={(v) => updateBox('auto_box', v)}
              description="Automatically determine box size based on protein structure"
            />
            {!settings.box.auto_box && (
              <>
                <Slider
                  label="Min Axis Size"
                  value={settings.box.min_axis}
                  min={10}
                  max={30}
                  onChange={(v) => updateBox('min_axis', v)}
                  description="Minimum box dimension in Angstroms"
                />
                <Slider
                  label="Max Axis Size"
                  value={settings.box.max_axis}
                  min={20}
                  max={50}
                  onChange={(v) => updateBox('max_axis', v)}
                  description="Maximum box dimension in Angstroms"
                />
              </>
            )}
          </div>
        </div>

        {/* Protein Cleaning */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-green-900/50 rounded-lg">
              <Droplets className="w-5 h-5 text-green-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">Protein Cleaning</h2>
              <p className="text-sm text-gray-400">Control what to keep during protein preparation</p>
            </div>
          </div>
          
          <div className="space-y-1">
            <Toggle
              label="Keep Waters"
              checked={settings.cleaning.keep_waters}
              onChange={(v) => updateCleaning('keep_waters', v)}
              description="Retain water molecules in the protein structure"
            />
            <Toggle
              label="Keep Ions"
              checked={settings.cleaning.keep_ions}
              onChange={(v) => updateCleaning('keep_ions', v)}
              description="Retain ion molecules (Na+, Cl-, etc.)"
            />
            <Toggle
              label="Keep Solvents"
              checked={settings.cleaning.keep_solvents}
              onChange={(v) => updateCleaning('keep_solvents', v)}
              description="Retain other solvent molecules"
            />
            <Toggle
              label="Keep Cofactors"
              checked={settings.cleaning.keep_cofactors}
              onChange={(v) => updateCleaning('keep_cofactors', v)}
              description="Retain essential cofactors (ATP, NAD, etc.)"
            />
          </div>
        </div>

        {/* Retention */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-yellow-900/50 rounded-lg">
              <Clock className="w-5 h-5 text-yellow-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">Data Retention</h2>
              <p className="text-sm text-gray-400">Control how long job data is kept</p>
            </div>
          </div>
          
          <RadioGroup
            label="Job Retention Policy"
            options={[
              { 
                value: 'permanent', 
                label: 'Keep permanently', 
                description: 'All job data is stored indefinitely' 
              },
              { 
                value: '7days', 
                label: 'Delete after 7 days', 
                description: 'Automatic cleanup of jobs older than 7 days' 
              },
              { 
                value: 'immediate', 
                label: 'Delete immediately after viewing', 
                description: 'Jobs are removed once results are viewed' 
              }
            ]}
            selected={settings.retention.retention}
            onChange={updateRetention}
          />
        </div>

        {/* Appearance */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-pink-900/50 rounded-lg">
              <Palette className="w-5 h-5 text-pink-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">Appearance</h2>
              <p className="text-sm text-gray-400">Customize the look and feel</p>
            </div>
          </div>
          
          <div className="bg-gray-750 rounded-lg p-4 text-center">
            <Palette className="w-10 h-10 text-gray-500 mx-auto mb-3" />
            <p className="text-gray-400 text-sm">Dark mode is currently the only theme.</p>
            <p className="text-gray-500 text-xs mt-1">Light mode support coming soon.</p>
          </div>
        </div>

        {/* About */}
        <div className="bg-gray-800 rounded-xl border border-gray-700 p-6">
          <div className="flex items-center gap-3 mb-6">
            <div className="p-2 bg-indigo-900/50 rounded-lg">
              <Info className="w-5 h-5 text-indigo-400" />
            </div>
            <div>
              <h2 className="text-lg font-semibold text-white">About</h2>
              <p className="text-sm text-gray-400">Application information</p>
            </div>
          </div>
          
          <div className="space-y-4">
            <div className="flex items-center justify-between py-2 border-b border-gray-700">
              <span className="text-gray-400">Version</span>
              <span className="text-white font-mono">0.1.0</span>
            </div>
            <div className="flex items-center justify-between py-2 border-b border-gray-700">
              <span className="text-gray-400">Framework</span>
              <span className="text-white">React + TypeScript</span>
            </div>
            <div className="flex items-center justify-between py-2 border-b border-gray-700">
              <span className="text-gray-400">Docking Engine</span>
              <span className="text-white">AutoDock Vina</span>
            </div>
            <div className="pt-2">
              <a
                href="https://github.com"
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center gap-2 text-blue-400 hover:text-blue-300 text-sm"
              >
                <ExternalLink className="w-4 h-4" />
                View on GitHub
              </a>
            </div>
          </div>
        </div>

        {/* Save Button (Bottom) */}
        <div className="flex justify-end pb-8">
          <button
            onClick={handleSave}
            disabled={saving}
            className="flex items-center gap-2 px-6 py-3 bg-blue-600 hover:bg-blue-700 text-white rounded-lg transition-colors disabled:opacity-50 font-medium"
          >
            {saving ? (
              <RefreshCw className="w-4 h-4 animate-spin" />
            ) : (
              <Save className="w-4 h-4" />
            )}
            {saving ? 'Saving...' : 'Save All Settings'}
          </button>
        </div>
      </div>
    </div>
  );
};

export default SettingsPage;
