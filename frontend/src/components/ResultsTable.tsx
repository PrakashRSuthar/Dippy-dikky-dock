// src/components/ResultsTable.tsx

interface Pose {
  pose: number;
  affinity: string;
  rmsd_lb: string;
  rmsd_ub: string;
}

interface ResultsTableProps {
  data: Pose[];
  selectedPose: number;
  onPoseSelect: (pose: number) => void;
}

export const ResultsTable = ({ data, selectedPose, onPoseSelect }: ResultsTableProps) => {
  return (
    <div className="overflow-x-auto">
      <table className="w-full text-sm">
        <thead>
          <tr className="border-b bg-gray-50">
            <th className="text-left py-3 px-4 font-medium text-gray-700">Pose</th>
            <th className="text-left py-3 px-4 font-medium text-gray-700">Affinity (kcal/mol)</th>
            <th className="text-left py-3 px-4 font-medium text-gray-700">RMSD l.b.</th>
            <th className="text-left py-3 px-4 font-medium text-gray-700">RMSD u.b.</th>
            <th className="text-left py-3 px-4 font-medium text-gray-700">Action</th>
          </tr>
        </thead>
        <tbody>
          {data.map((pose, idx) => (
            <tr 
              key={idx}
              className={`border-b hover:bg-gray-50 ${
                selectedPose === idx ? 'bg-blue-50 border-blue-200' : ''
              }`}
            >
              <td className="py-3 px-4 font-medium">
                {selectedPose === idx && <span className="text-blue-600">► </span>}
                {pose.pose}
              </td>
              <td className="py-3 px-4">
                <span className={`font-medium ${
                  parseFloat(pose.affinity) <= -8 ? 'text-green-600' :
                  parseFloat(pose.affinity) <= -6 ? 'text-yellow-600' : 'text-red-600'
                }`}>
                  {pose.affinity}
                </span>
              </td>
              <td className="py-3 px-4 text-gray-600">{pose.rmsd_lb}</td>
              <td className="py-3 px-4 text-gray-600">{pose.rmsd_ub}</td>
              <td className="py-3 px-4">
                <button
                  onClick={() => onPoseSelect(idx)}
                  className={`px-3 py-1 rounded-md text-xs font-medium ${
                    selectedPose === idx
                      ? 'bg-blue-600 text-white'
                      : 'bg-gray-100 text-gray-700 hover:bg-gray-200'
                  }`}
                >
                  {selectedPose === idx ? 'Selected' : 'View'}
                </button>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};
