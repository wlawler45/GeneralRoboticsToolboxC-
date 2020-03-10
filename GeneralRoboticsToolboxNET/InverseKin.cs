// Copyright 2020 Rensselaer Polytechnic Institute
//                Wason Technology, LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace GeneralRoboticsToolbox
{
    public class NormalizeJoints
    {
        public NormalizeJoints() { }

        public Robot Robot { get; set; }
        public double[] Last_joints { get; set; }
        private bool check_limits { get; set; }
        private bool use_last_joints { get; set; }
        private double current_compare { get; set; }
        private double[] current_compare2 { get; set; }


        public NormalizeJoints(Robot robot, double[] last_joints = default(double[]))
        {
            Robot = robot;
            Last_joints = last_joints;
            check_limits = Robot.Joint_upper_limit != null && Robot.Joint_lower_limit != null;
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
        public double[] FindNormalizedJoints(int joint_index, double[] thetas)
        {
            List<double> theta_normed = new List<double>();
            foreach (double t1 in thetas)
            {
                double t3 = Normalize(joint_index, t1);
                if (t3 != default(double)) theta_normed.Add(t3);
            }

            if (!use_last_joints || theta_normed.Count < 2)
            {
                double[] output = new double[theta_normed.Count];
                for (int i = 0; i < theta_normed.Count; i++)
                {
                    output[i] = theta_normed[i];
                }
                return output;
            }
            // theta_last = Vector<double>.Build.Dense(joint_index.Length);
            Vector<double> theta_normed_vec = Vector<double>.Build.DenseOfArray(theta_normed.ToArray());

            double theta_last = Last_joints[joint_index];

            Vector<double> theta_dist_vec = theta_normed_vec.Subtract(theta_last);
            theta_dist_vec = theta_dist_vec.PointwiseAbs();
            //double[] theta_dist = theta_dist_vec.ToArray();
            double[] theta_dist = theta_dist_vec.ToArray();
            List<double> theta_ret1 = new List<double>();
            foreach (double theta in theta_normed)
            {
                current_compare = theta;
                if (Math.Abs(current_compare - theta_last) < Math.PI / 2.0) theta_ret1.Add(theta);
            }
            if (theta_ret1.Count == 1)
            {
                return theta_ret1.ToArray();
            }
            if (theta_ret1.Count == 0)
            {
                theta_ret1 = theta_normed;
            }
            double[] returned = new double[theta_dist.Length];
            int index = 0;
            List<double> theta_dist_list = theta_dist.ToList();
            for (int z = 0; z < theta_dist.Length; z++)
            {
                index = theta_dist_list.IndexOf(theta_dist_list.Min());
                returned[z] = theta_normed[index];
                theta_dist_list.RemoveAt(index);
            }
            return returned;


        }

        public double[][] FindNormalizedJoints(int[] joint_index, double[] thetas)
        {
            List<List<double>> theta_normed = new List<List<double>>();

            foreach (double t1 in thetas)
            {
                List<double> theta_normed_intermed = new List<double>();
                for (int i = 0; i < joint_index.Length; i++)
                {
                    double t4 = Normalize(joint_index[i], t1);
                    if (t4 != default(double)) theta_normed_intermed.Add(t4);
                }
                theta_normed.Add(theta_normed_intermed);
            }



            if (!use_last_joints || theta_normed.Count < 2)
            {

                double[][] output = theta_normed.Select(a => a.ToArray()).ToArray();
                return output;
            }
            Vector<double> theta_last = Vector<double>.Build.Dense(joint_index.Length);

            for (int i = 0; i < joint_index.Length; i++)
            {

                theta_last[i] = Last_joints[joint_index[i]];
            }

            double[] theta_dist_arr = new double[joint_index.Length];

            List<List<double>> theta_ret = new List<List<double>>();
            for (int i = 0; i < joint_index.Length; i++)
            {

                Vector<double> theta_normed_vec = Vector<double>.Build.DenseOfArray(theta_normed[i].ToArray());
                Vector<double> theta_dist_vec = theta_normed_vec.Subtract(theta_last[i]);
                theta_dist_arr[i] = theta_dist_vec.L2Norm();

                bool passed = true;
                foreach (double t in theta_dist_vec)
                {
                    if (!(Math.Abs(t) < Math.PI / 2.0)) passed = false;
                }
                if (passed) theta_ret.Add(theta_normed[i]);
            }

            /*if (joint_index.Length== 1)
            {
                
            }
            else
            {*/

            //Vector<double> theta_dist = theta_dist_vec.PointwiseAbs();
            //just need array of pointwise magnitude values 
            // double[] theta_dist = theta_dist_vec.ToArray();
            //}


            /*foreach(double theta in theta_normed)
            {
                current_compare = theta; 
                if (theta_last.ForAll(theta_vec_test)) theta_ret1.Add(theta);
            }*/
            if (theta_ret.Count == 1)
            {
                return theta_ret.Select(a => a.ToArray()).ToArray();
            }
            if (theta_ret.Count == 0)
            {
                theta_ret = theta_normed;
            }
            //double[] returned = new double[theta_dist.Count];
            int index = 0;
            List<double> theta_dist_list = theta_dist_arr.ToList();
            List<List<double>> returned = new List<List<double>>();
            for (int z = 0; z < theta_dist_list.Count; z++)
            {
                index = theta_dist_list.IndexOf(theta_dist_list.Min());
                returned[z] = theta_normed[index];
                theta_dist_list.RemoveAt(index);
            }
            return returned.Select(a => a.ToArray()).ToArray();

        }
    }
    public class InverseKin
    {

    
        public static double[][] robot6_sphericalwrist_invkin(Robot robot, Transform desired_pose, double[] last_joints = default(double[]))
        {

            if (robot.R_tool != default(Matrix<double>) && robot.P_tool != default(Vector<double>)){
                
                Matrix<double> transposed= robot.R_tool.Transpose();
                desired_pose.R = desired_pose.R.Multiply(transposed);
                desired_pose.P = desired_pose.P-desired_pose.R*robot.P_tool;
            }
            Matrix<double> H = robot.H.Clone();
            Matrix<double> P= robot.P.Clone();
            List<double[]> theta_v = new List<double[]>();
            Vector<double> zeros = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 0.0 });
            Vector<double> ex = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0.0, 0.0 });
            Vector<double> ey = Vector<double>.Build.DenseOfArray(new[] { 0.0, 1.0, 0.0 });
            Vector<double> ez = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 1.0 });
            if (!P.Column(4).Equals(zeros))
            {
                double P4_d = P.Column(4) * H.Column(3);
                if (!(P.Column(4).Subtract(P4_d * H.Column(3)).Equals(zeros)) )throw new ArgumentException(String.Format("Robot may not have spherical wrist"));
                P.SetColumn(3, P.Column(3) + P.Column(4));
                P.SetColumn(4, zeros);
            }
            if (!P.Column(5).Equals(zeros))
            {
                double P5_d = P.Column(5) * H.Column(5);
                if (!(P.Column(5).Subtract(P5_d * H.Column(5)).Equals(zeros))) throw new ArgumentException(String.Format("Robot may not have spherical wrist"));
                P.SetColumn(6, P.Column(6) + P.Column(5));
                P.SetColumn(5, zeros);
            }

            double d1 = ey*(P.Column(1).Add(P.Column(2)).Add(P.Column(3)));
            Vector<double> v1 = desired_pose.P.Subtract( desired_pose.R * P.Column(6));
            double[] Q1 = Functions.Subproblem4(ey, v1, -H.Column(0), d1);
            NormalizeJoints normalize = new NormalizeJoints(robot, last_joints);
            
            double[] first_normalize = normalize.FindNormalizedJoints(0, Q1);
            foreach (double q1 in first_normalize)
            {
                Matrix<double> R01 = Functions.Rot(H.Column(0), q1);
                Vector<double> p26_q1 = R01.TransposeThisAndMultiply(desired_pose.P - desired_pose.R * P.Column(6)).Subtract(P.Column(0).Add(P.Column(1)));
                double d3 = p26_q1.L2Norm();
                Vector<double> v3 = P.Column(2);
                Vector<double> p3 = P.Column(3);
                double[] Q3 = Functions.Subproblem3(p3, v3, H.Column(2), d3);
                double[] second_normalize = normalize.FindNormalizedJoints(2, Q3);
                foreach (double q3 in second_normalize)
                {
                    Matrix<double> R23 = Functions.Rot(H.Column(2), q3);
                    Vector<double> v2 = p26_q1;
                    Vector<double> p2 = P.Column(2) + R23 * P.Column(3);
                    double q2 = Functions.Subproblem1(p2, v2, H.Column(1));
                    double[] q2_array = normalize.FindNormalizedJoints(1, new[] { q2 });
                    if (q2_array.Length == 0) continue;
                    q2 = q2_array[0];
                    Matrix<double> R12 = Functions.Rot(H.Column(1), q2);
                    Matrix<double> R03 = R01 * R12 * R23;
                    Matrix<double> R36 = R03.Transpose() * desired_pose.R;
                    Vector<double> v4 = R36 * H.Column(5);
                    Vector<double> p4 = H.Column(5);
                    double[] Q4_Q5 = Functions.Subproblem2(p4, v4, H.Column(3), H.Column(4));
                    double[][] third_normalize = normalize.FindNormalizedJoints(new[] { 3,4}, Q4_Q5);
                    Console.WriteLine("size of qs {0}, {1}", third_normalize[0].Length, third_normalize[1].Length);
                    int minoftwo = Math.Min(third_normalize[0].Length, third_normalize[1].Length);
                    Console.WriteLine(minoftwo);
                    for (int q=0; q < minoftwo; q++)
                    {
                        Matrix<double> R35 = Functions.Rot(H.Column(3), third_normalize[0][q])*Functions.Rot(H.Column(4), third_normalize[1][q]);
                        Matrix<double> R05 = R03 * R35;
                        Matrix<double> R56 = R05.Transpose() * desired_pose.R;
                        Vector<double> v6 = R56 * H.Column(4);
                        Vector<double> p6 = H.Column(4);
                        double q6 = Functions.Subproblem1(p6, v6, H.Column(5));
                        double[] q6_array = normalize.FindNormalizedJoints(5, new[] { q6 });
                        if (q6_array.Length == 0) continue;
                        q6 = q6_array[0];
                        
                        double[] theta_v_entry = new double[] { q1, q2, q3, third_normalize[0][q], third_normalize[1][q], q6 };
                        theta_v.Add(theta_v_entry);


                    }


                }
            }

            if (last_joints != default(double[]))
            {
                double[] theta_dist_arr = new double[theta_v.Count];


                for (int i = 0; i < theta_v.Count; i++)
                {

                    Vector<double> theta_v_vec = Vector<double>.Build.DenseOfArray(theta_v[i]);
                    Vector<double> theta_dist_vec = theta_v_vec.Subtract(last_joints[i]);
                    theta_dist_arr[i] = theta_dist_vec.L2Norm();
                }
                int index = 0;
                List<double> theta_dist_list = theta_dist_arr.ToList();
                List<double[]> returned = new List<double[]>();
                for (int z = 0; z < theta_dist_list.Count; z++)
                {
                    index = theta_dist_list.IndexOf(theta_dist_list.Min());
                    returned[z] = theta_v[index];
                    theta_dist_list.RemoveAt(index);
                }
                return returned.ToArray();
            }
            else
            {
                return theta_v.ToArray();
            }
        }

    }
    
}
