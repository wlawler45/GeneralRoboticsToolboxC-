using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace TestGeneralRoboticsToolboxNET
{
    public class InverseKin
    {
        public InverseKin() { }
        public Matrix<double> R { get; set; }
        public Vector<double> P { get; set; }
        public string Parent_frame_id { get; set; }
        public string Child_frame_id { get; set; }

        public void normalize_joints(Robot robot, Vector<double> last_joints)
        {
            //if (r.ColumnCount != 3 || r.RowCount != 3) throw new ArgumentException(String.Format("Rotation Matrix for transform not Acceptable"));
            //if (p.Count != 3) throw new ArgumentException(String.Format("Position Vector for transform not Acceptable"));
            //R = r;
            //P = p;
            //Parent_frame_id = parent_frame_id;
            //Child_frame_id = child_frame_id;
        }
    }
}
