// src/components/ResultsChart.tsx

interface Pose {
  pose: number;
  affinity: string;
  rmsd_lb: string;
  rmsd_ub: string;
}

interface ResultsChartProps {
  data: Pose[];
  selectedPose: number;
}

export const ResultsChart = ({ data, selectedPose }: ResultsChartProps) => {
  const maxAffinity = Math.max(...data.map(p => Math.abs(parseFloat(p.affinity))));
  
  return (
    <div className="w-full h-full flex items-end space-x-1 px-4 py-4">
      {data.slice(0, 20).map((pose, idx) => {
        const height = (Math.abs(parseFloat(pose.affinity)) / maxAffinity) * 100;
        const isSelected = selectedPose === idx;
        
        return (
          <div
            key={idx}
            className="flex flex-col items-center flex-1 min-w-0"
          >
            <div className="text-xs text-gray-600 mb-1 transform -rotate-45 origin-bottom-left">
              {pose.affinity}
            </div>
            <div
              className={`w-full rounded-t transition-all duration-200 ${
                isSelected 
                  ? 'bg-blue-600' 
                  : parseFloat(pose.affinity) <= -8 
                    ? 'bg-green-500' 
                    : parseFloat(pose.affinity) <= -6 
                      ? 'bg-yellow-500' 
                      : 'bg-red-500'
              }`}
              style={{ height: `${height}%`, minHeight: '4px' }}
            />
            <div className="text-xs text-gray-500 mt-1">{pose.pose}</div>
          </div>
        );
      })}
    </div>
  );
};
