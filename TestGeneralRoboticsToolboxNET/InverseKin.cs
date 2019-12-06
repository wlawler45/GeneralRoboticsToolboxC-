using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace TestGeneralRoboticsToolboxNET
{
    public class NormalizeJoints
    {
        public NormalizeJoints() { }
        
        public Robot Robot { get; set; }
        public double[] Last_joints { get; set; }
        private bool check_limits { get; set; }
        private bool use_last_joints { get; set; }
        


        public NormalizeJoints(Robot robot, double[] last_joints=default(double[]))
        {
            Robot = robot;
            Last_joints = last_joints;
            check_limits = Robot.Joint_upper_limit != null && Robot.Joint_lower_limit!=null;
            use_last_joints = last_joints != default(double[]);
            //if (r.ColumnCount != 3 || r.RowCount != 3) throw new ArgumentException(String.Format("Rotation Matrix for transform not Acceptable"));
            //if (p.Count != 3) throw new ArgumentException(String.Format("Position Vector for transform not Acceptable"));
            //R = r;
            //P = p;
            //Parent_frame_id = parent_frame_id;
            //Child_frame_id = child_frame_id;
        }
        public double Normalize(int joint_index, double theta)
        {
            double u = Robot.Joint_upper_limit[joint_index];
            double l = Robot.Joint_lower_limit[joint_index];
            if (check_limits)
            {

                if (!(l < theta && theta < u))
                {
                    Vector<double> a = Vector<double>.Build.DenseOfArray(new[] { -1.0, 1.0 });
                    a = a * 2 * Math.PI;
                    Vector<double> b = a.Add(theta);
                    int index = -1;
                    for (int i = 0; i < b.Count; i++)
                    {
                        if (l < b[i] && b[i] < u)
                        {
                            index = i;
                            break;
                        }
                    }
                    if (index == -1)
                    {
                        return default(double);
                    }
                    else
                    {
                        theta = theta + a[index];
                    }
                }
            }
            if (use_last_joints)
            {
                double diff = Last_joints[joint_index] - theta;

                int n_diff = (int)Math.Floor(diff / (Math.PI * 2.0));
                double r_diff = diff % (Math.PI * 2.0);
                if (r_diff > Math.PI)
                {
                    n_diff += 1;
                }
                if (Math.Abs(n_diff) > 0)
                {
                    if (!check_limits)
                    {
                        theta += 2 * Math.PI * n_diff;
                    }
                    else
                    {
                        double theta_vs = theta + 2 * Math.PI;

                        Vector<double> theta_v = Vector<double>.Build.Dense(Math.Abs(n_diff));
                        for (int i = 0; i < theta_v.Count; i++)
                        {
                            if (n_diff > 0)
                            {
                                theta_v[i] = n_diff - i;
                            }
                            else if (n_diff < 0) {
                                theta_v[i] = n_diff + i;
                            }
                        }
                        theta_v = theta_v * 2 * Math.PI;
                        theta_v = theta_v.Add(theta);
                        int index = -1;
                        for (int i = 0; i < theta_v.Count; i++)
                        {
                            if (l < theta_v[i] && theta_v[i] < u)
                            {
                                index = i;
                                break;
                            }
                        }
                        if (index != -1)
                        {
                            theta = theta_v[index];
                        }
                    }
                }


            }
            return theta;
        }
        public double[] FindNormalizedJoints(int[] joint_index, double[] thetas)
        {
            List<double> theta_normed = new List<double>();
            if (joint_index.Length == 1)
            {
                foreach(double t1 in thetas)
                {
                    double t3 = Normalize(joint_index[0], t1);
                    if (t3 != default(double)) theta_normed.Add(t3);
                }
            }
            else
            {
                
                for(int i = 0; i < joint_index.Length; i++)
                {
                    double t4=Normalize(joint_index[i], thetas[i]);
                    if (t4 != default(double)) theta_normed.Add(t4);
                }
                
                
                
            }


        }
    }
}
