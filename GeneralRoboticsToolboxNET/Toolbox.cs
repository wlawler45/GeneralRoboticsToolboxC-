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
using System.Diagnostics;
using System.Text;
using System.Numerics;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Spatial.Euclidean;

namespace GeneralRoboticsToolbox
{
    public static class Functions
    {        
        public static Matrix<double> Hat(Vector<double> k)
        {
            Matrix<double> khat = Matrix<double>.Build.Dense(3, 3);
            khat[0, 1] = -k[2];
            khat[0, 2] = k[1];
            khat[1, 0] = k[2];
            khat[1, 2] = -k[0];
            khat[2, 0] = -k[1];
            khat[2, 1] = k[0];
            return khat;
        }

        public static Vector<double> Invhat(Matrix<double> khat)
        {
            double[] inv = { (-khat[1, 2] + khat[2, 1]), (khat[0, 2] - khat[2, 0]), (-khat[0, 1] + khat[1, 0]) };
            Vector<double> output = Vector<double>.Build.DenseOfArray(inv);
            output /= 2;
            return output;
        }
        public static Matrix<double> Rot(Vector<double> k, double theta)
        {
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(3);
            Matrix<double> khat = Hat(k);
            Matrix<double> khat2 = khat.Multiply(khat);
            return I + Math.Sin(theta) * khat + (1.0 - Math.Cos(theta)) * khat2;
        }
        public static Tuple<Vector<double>, double> R2rot(Matrix<double> R)
        {
            Matrix<double> R1 = R - R.Transpose();
            //double sin_theta = R1.L2Norm() / Math.Sqrt(8);
            double sin_theta = R1.L2Norm() / 2.0;
            double cos_theta = (R.Trace() - 1.0) / 2.0;
            double theta = Math.Atan2(sin_theta, cos_theta);
            Vector<double> k;
            if (sin_theta < (1e-6))
            {
                if (cos_theta > 0)
                {
                    k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 1.0 });

                    return new Tuple<Vector<double>, double>(k, 0);
                }
                else
                {
                    Matrix<double> eye = Matrix<Double>.Build.DenseIdentity(3);
                    Matrix<double> B = (1.0 / 2.0) * (R + eye);
                    k = Vector<double>.Build.DenseOfArray(new[] { Math.Sqrt(B[0, 0]), Math.Sqrt(B[1, 1]), Math.Sqrt(B[2, 2]) });
                    if (Math.Abs(k[0]) > (1e-6))
                    {
                        k[1] = k[1] * Math.Sign(B[0, 1] / k[0]);
                        k[2] = k[2] * Math.Sign(B[0, 2] / k[0]);
                    }
                    else if (Math.Abs(k[1]) > (1e-6))
                    {
                        k[2] = k[2] * Math.Sign(B[0, 2] / k[1]);
                    }
                    return new Tuple<Vector<double>, double>(k, Math.PI);

                }
            }
            k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 0.0 });
            Vector<double> inv = Invhat(R1);

            for (int i = 0; i < k.Count; i++)
            {
                k[i] = inv[i] / (2.0 * sin_theta);
            }

            return new Tuple<Vector<double>, double>(k, theta);
        }
        public static Matrix<double> Screw_matrix(Vector<double> r)
        {
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(6);
            Matrix<double> hat = Hat(r);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 3; j < 6; j++)
                {
                    I[i, j] = hat[i, j - 3];
                }
            }
            return I;
        }
        public static Matrix<double> Q2R(Vector<double> q)
        {
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(3);

            Matrix<double> hat = Hat(q.SubVector(1, 3));
            
            Matrix<double> qhat2 = hat.Multiply(hat);
            return I + 2 * q[0] * hat + 2 * qhat2;

        }
        public static Vector<double> R2Q(Matrix<double> R)
        {
            double tr = R.Trace();
            Vector<double> q;
            if (tr > 0)
            {
                double S = 2 * Math.Sqrt(tr + 1);

                q = Vector<double>.Build.DenseOfArray(new[] {(0.25 * S),
                      ((R[2, 1] - R[1, 2]) / S),
                      ((R[0, 2] - R[2, 0]) / S),
                      ((R[1, 0] - R[0, 1]) / S)});

            }
            else if (R[0, 0] > R[1, 1] && R[0, 0] > R[2, 2])
            {
                double S = 2 * Math.Sqrt(1 + R[0, 0] - R[1, 1] - R[2, 2]);
                q = Vector<double>.Build.DenseOfArray(new[] {((R[2, 1] - R[1, 2]) / S),
                      (0.25 * S),
                      ((R[0, 1] + R[1, 0]) / S),
                      ((R[0, 2] + R[2, 0]) / S)});
            }
            else if (R[1, 1] > R[2, 2])
            {
                double S = 2 * Math.Sqrt(1 - R[0, 0] + R[1, 1] - R[2, 2]);
                q = Vector<double>.Build.DenseOfArray(new[] {((R[0, 2] - R[2, 0]) / S),
                      ((R[0, 1] + R[1, 0]) / S),
                      (0.25 * S),
                      ((R[1, 2] + R[2, 1]) / S) });

            }
            else
            {
                double S = 2 * Math.Sqrt(1 - R[0, 0] - R[1, 1] + R[2, 2]);
                q = Vector<double>.Build.DenseOfArray(new[] {((R[1, 0] - R[0, 1]) / S),
                      ((R[0, 2] + R[2, 0]) / S),
                      ((R[1, 2] + R[2, 1]) / S),
                      (0.25 * S) });
            }
            return q;

        }
        public static Tuple<Vector<double>, double> Q2Rot(Vector<double> q)
        {
            double theta = 2 * Math.Acos(q[0]);
            Console.WriteLine(theta);
            Vector<double> k;
            if (Math.Abs(theta) < 1e-6)
            {
                k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 1.0 });

                return new Tuple<Vector<double>, double>(k, 0);
            }
            k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 1.0 });
            for (int i = 1; i < k.Count + 1; i++)
            {
                k[i - 1] = q[i] / Math.Sin(theta / 2.0);
            }
            return new Tuple<Vector<double>, double>(k, theta);
        }
        public static Vector<double> Rot2Q(Vector<double> k, double theta)
        {
            double s1 = Math.Sin(theta / 2.0);
            double c1 = Math.Cos(theta / 2.0);            
            Vector<double> output = Vector<double>.Build.DenseOfArray(new[] { c1, k[0]*s1, k[1]*s1, k[2]*s1 });
            return output;
        }
        public static Vector<double> Quatcomplement(Vector<double> q)
        {
            Vector<double> output = Vector<double>.Build.DenseOfArray(new[] { q[0], -1 * q[1], -1 * q[2], -1 * q[3] });
            return output;
        }
        public static Matrix<double> Quatproduct(Vector<double> q)
        {
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(3);
            Matrix<double> Q = Matrix<double>.Build.Dense(4, 4);
            Matrix<double> hats = Hat(q.SubVector(1, 3));
            Q[0, 0] = q[0];

            for (int i = 1; i < 4; i++)
            {
                Q[0, i] = -q[i];
            }
            for (int i = 1; i < 4; i++)
            {
                Q[i, 0] = q[i];
            }
            for (int i = 1; i < 4; i++)
            {
                for (int j = 1; j < 4; j++)
                {
                    Q[i, j] = q[0] * I[i - 1, j - 1] + hats[i - 1, j - 1];
                }

            }
            return Q;
        }
        public static Matrix<double> Quatjacobian(Vector<double> q)
        {
            Matrix<double> I = Matrix<double>.Build.DenseIdentity(3);
            Matrix<double> J = Matrix<double>.Build.Dense(4, 3);
            Matrix<double> hats = Hat(q.SubVector(1, 3));
            for (int i = 0; i < 3; i++)
            {
                J[0, i] = 0.5 * -q[i + 1];
            }
            for (int i = 1; i < 4; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    J[i, j] = 0.5 * (q[0] * I[i - 1, j] - hats[i - 1, j]);
                }

            }
            return J;
        }
        public static Matrix<double> rpy2R(Vector<double> rpy)
        {
            Vector<double> k;
            Matrix<double> rotation1 = Rot(k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0.0, 1.0 }), rpy[2]);
            Matrix<double> rotation2 = Rot(k = Vector<double>.Build.DenseOfArray(new[] { 0.0, 1.0, 0.0 }), rpy[1]);
            Matrix<double> rotation3 = Rot(k = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0.0, 0.0 }), rpy[0]);
            return rotation1.Multiply(rotation2).Multiply(rotation3);
        }
        public static Vector<double> R2rpy(Matrix<double> R)
        {
            //if(R.Column(0).L2Norm())
            Vector<double> output;
            double r = Math.Atan2(R[2, 1], R[2, 2]);
            double y = Math.Atan2(R[1, 0], R[0, 0]);
            //Console.WriteLine(R);
            double normie = (R.SubMatrix(2, 1, 1, 2)).L2Norm();
            double p = Math.Atan2(-R[2, 0], normie);
            return output = Vector<double>.Build.DenseOfArray(new[] { r, p, y });
        }

        public static Transform Fwdkin(Robot robot, double[] theta)
        {
            
            if(robot.Joint_lower_limit!=null && robot.Joint_upper_limit != null)
            {
                for(int i = 0; i < theta.Length; i++)
                {
                    if(!(robot.Joint_lower_limit[i]<theta[i] && theta[i]<robot.Joint_upper_limit[i])) throw new ArgumentException(String.Format("Joint angle for joint {0} out of range",i));
                }
            }
            Vector<double> p;
            p = robot.P.Column(0);
            Matrix<double> R = Matrix<double>.Build.DenseIdentity(3);
            for(int i = 0; i < robot.Joint_type.Length; i++)
            {
                if(robot.Joint_type[i]== JointType.Revolute || robot.Joint_type[i] == JointType.MobileOrientation)
                {
                    R = R * Rot(robot.H.Column(i), theta[i]);
                }else if(robot.Joint_type[i] == JointType.Prismatic || robot.Joint_type[i] == JointType.MobileTranslation)
                {
                    p = p + theta[i] * R * robot.H.Column(i);
                }
                p = p + R * robot.P.Column(i + 1);
            }
            
            //Reshape here not added
            if (robot.R_tool != null && robot.P_tool != null)
            {
                
                p = p + R.Multiply(robot.P_tool);
                R = R*robot.R_tool;
                
                //p = p + R*robot.P_tool;
                //R = R*robot.R_tool;

            }
            Transform output = new Transform(R, p);
            
            return output;

        }

        public static Matrix<double> Robotjacobian(Robot robot, double[] theta)
        {

            if (robot.Joint_lower_limit != null && robot.Joint_upper_limit != null)
            {
                for (int k = 0; k < theta.Length; k++)
                {
                    if (!(robot.Joint_lower_limit[k] < theta[k] && theta[k] < robot.Joint_upper_limit[k])) throw new ArgumentException(String.Format("Joint angles out of range"));
                }
            }

            Matrix<double> hi = Matrix<double>.Build.Dense(robot.H.RowCount, robot.H.ColumnCount);
            Matrix<double> p0i = Matrix<double>.Build.Dense(robot.P.RowCount, robot.P.ColumnCount);
            Vector<double> p;
            p = robot.P.Column(0);
            Matrix<double> R = Matrix<double>.Build.DenseIdentity(3);

            p0i.SetColumn(0,p);
            


            for (int k = 0; k < robot.Joint_type.Length; k++)
            {
                if (robot.Joint_type[k] == 0 || robot.Joint_type[k] == JointType.MobileOrientation)
                {
                    R = R * Rot(robot.H.Column(k), theta[k]);
                }
                else if (robot.Joint_type[k] == JointType.Prismatic || robot.Joint_type[k] == JointType.MobileTranslation)
                {
                    p = p + theta[k] * R * robot.H.Column(k);
                }
                p = p + R * robot.P.Column(k + 1);
                p0i.SetColumn(k + 1, p);
                hi.SetColumn(k, R * robot.H.Column(k));

            }
            Vector<double> p0T = p0i.Column(robot.Joint_type.Length);
            

            if (robot.P_tool != null)
            {
                p0T = p0T + R * robot.P_tool;
            }
            Matrix<double> J= Matrix<double>.Build.Dense(6, robot.Joint_type.Length);
            int i = 0;
            int j = 0;
            while (i < robot.Joint_type.Length)
            {
                if (robot.Joint_type[i] == JointType.Revolute)
                {
                    J.SetColumn(j, 0, 3, hi.Column(i));
                    J.SetColumn(j, 3, 3, Hat(hi.Column(i)) * (p0T - p0i.Column(i)));
                }else if(robot.Joint_type[i] == JointType.Prismatic)
                {
                    J.SetColumn(j, 3, 3, hi.Column(i));
                }else if(robot.Joint_type[i] == JointType.MobileTranslation)
                {
                    J.SetColumn(j, 3, 3, Rot(hi.Column(i + 2), theta[i + 2]) * hi.Column(i));
                    J.SetColumn(j + 1, 0, 3, hi.Column(i + 2));
                    J.SetColumn(j+1, 3, 3, Hat(hi.Column(i+2)) * (p0T - p0i.Column(i+2)));
                    J = J.SubMatrix(0, J.RowCount, 0, J.ColumnCount - 1);
                    i+=2;
                    j++;
                }
                i++;
                j++;
            }
            return J;


            

        }
        public static Vector<double> Cross(Vector<double> left, Vector<double> right)
        {
            if ((left.Count != 3 || right.Count != 3))
            {
                string message = "Vectors must have a length of 3.";
                throw new Exception(message);
            }
            Vector<double> result =  Vector<double>.Build.Dense(3);
            result[0] = left[1] * right[2] - left[2] * right[1];
            result[1] = -left[0] * right[2] + left[2] * right[0];
            result[2] = left[0] * right[1] - left[1] * right[0];

            return result;
        }
        public static double Subproblem0(Vector<double> p, Vector<double> q, Vector<double> k)
        {
            double min = double.Epsilon;
            if (!(k * p < min && k * q < min)) throw new ArgumentException(String.Format("k must be perpendicular to p and q"));
            Vector<double> ep = p / p.L2Norm();
            Vector<double> eq = q / q.L2Norm();
            double theta = 2 * Math.Atan2((ep - eq).L2Norm(), (ep + eq).L2Norm());
            if (k * Cross(p, q) < 0) return -theta;
            return theta;
        }

        public static double Subproblem1(Vector<double> p, Vector<double> q, Vector<double> k)
        {
            double min = double.Epsilon;
            Vector<double> minus = p - q;
            if (minus.L2Norm() < Math.Sqrt(min)) return 0.0;
            k = k / k.L2Norm();
            Vector<double> pp = p - (p * k * k);
            Vector<double> qp = q - (q * k * k);

            Vector<double> epp = pp / pp.L2Norm();
            Vector<double> eqp = qp / qp.L2Norm();

            double theta = Subproblem0(epp, eqp, k);
            if (Math.Abs(p.L2Norm() - q.L2Norm()) > (p.L2Norm() * (1e-2))) Console.WriteLine("WARNING:||p|| and ||q|| must be the same!!!");
            
            return theta;
        }
        public static double[] Subproblem2(Vector<double> p, Vector<double> q, Vector<double> k1, Vector<double> k2)
        {
            double min = double.Epsilon;
            double k12 = k1.DotProduct(k2);
            double pk = p.DotProduct(k2);
            double qk = q.DotProduct(k1);
            if (Math.Abs(1 - Math.Pow(k12, 2)) < min)
            {
                Console.WriteLine("WARNING:No solution found k1 != k2");
                return new double[0];
            }
            Matrix<double> amatrix = Matrix<double>.Build.Dense(2, 2);
            amatrix[0, 0] = k12;
            amatrix[0, 1] = -1;
            amatrix[1, 0] = -1;
            amatrix[1, 1] = k12;
            Vector<double> avector = Vector<double>.Build.DenseOfArray(new[] { pk, qk });
            Vector<double> a = amatrix.Multiply(avector / (Math.Pow(k12, 2) - 1));
            double bb = p.DotProduct(p) - a.DotProduct(a) - 2 * a[0] * a[1] * k12;
            if (Math.Abs(bb) < min) bb = 0;
            if (bb < 0)
            {
                Console.WriteLine("WARNING:No solution found no intersection found between cones");
                return new double[0];
            }
            double gamma = Math.Sqrt(bb) / Cross(k1, k2).L2Norm();
            if (Math.Abs(gamma) < min)
            {
                Matrix<double> cmalt = Matrix<double>.Build.DenseOfRowVectors(k1, k2, Cross(k1, k2));
                cmalt = cmalt.Transpose();
                Vector<double> c1vecalt = Vector<double>.Build.Dense(3);
                c1vecalt[0] = a[0];
                c1vecalt[1] = a[1];
                c1vecalt[2] = gamma;

                Vector<double> c1alt = cmalt.Multiply(c1vecalt);
                double theta2 = Subproblem1(k2, p, c1alt);
                double theta1 = -Subproblem1(k1, q, c1alt);
                double[] thetasfirst = { theta1, theta2 };
               
                
                return thetasfirst;
            }
            Matrix<double> cm = Matrix<double>.Build.DenseOfRowVectors(k1, k2, Cross(k1, k2));
            cm = cm.Transpose();
            Vector<double> c1vec = Vector<double>.Build.Dense(3);
            c1vec[0] = a[0];
            c1vec[1] = a[1];
            c1vec[2] = gamma;
            Vector<double> c2vec = Vector<double>.Build.Dense(3);
            c2vec[0] = a[0];
            c2vec[1] = a[1];
            c2vec[2] = -gamma;
            Vector<double> c1 = cm.Multiply(c1vec);
            Vector<double> c2 = cm.Multiply(c2vec);

            double theta1_1 = -Subproblem1(q, c1, k1);
            double theta1_2 = -Subproblem1(q, c2, k1);
            double theta2_1 = Subproblem1(p, c1, k2);
            double theta2_2 = Subproblem1(p, c2, k2);
            double[] thetas = new double[4] { theta1_1, theta2_1, theta1_2, theta2_2 };
            return thetas;
            //NOTE:THIS DOES NOT RETURN MULTIDIM ARRAY LIKE PYTHON VERSION
        }
        public static double[] Subproblem3(Vector<double> p, Vector<double> q, Vector<double> k, double d)
        {

            Vector<double> pp = p - ((p * k) * k);
            Vector<double> qp = q - ((q * k) * k);
            double dpsq = Math.Pow(d, 2) - Math.Pow((k.DotProduct(p + q)), 2);
            double bb = -(pp.DotProduct(pp) + qp.DotProduct(qp) - dpsq) / (2 * pp.L2Norm() * qp.L2Norm());
            if (dpsq < 0 || Math.Abs(bb) > 1)
            {
                Console.WriteLine("No solution no rotation can achieve specified distance");
                return new double[0];
            }
            double theta = Subproblem1(pp / pp.L2Norm(), qp / qp.L2Norm(), k);
            double phi = Math.Acos(bb);
            Console.WriteLine("theta, phi {0}, {1}", theta, phi);
            if (Math.Abs(phi) > 0)
            {
                return new double[2] { theta + phi, theta - phi };
            }
            else
            {
                return new double[1] { theta };
            }

        }
        public static double[] Subproblem4(Vector<double> p, Vector<double> q, Vector<double> k, double d)
        {

            Matrix<double> hatted = Hat(k);
            double a = hatted.LeftMultiply(p).DotProduct(q);

            double b = -1*hatted.LeftMultiply(hatted.LeftMultiply(q)).DotProduct(p);
            double c = d - (p.DotProduct(q) - b);
            double phi = Math.Atan2(b, a);
            Vector<double> dprep = Vector<double>.Build.Dense(2);
            dprep[0] = a;
            dprep[1] = b;
            d = c / dprep.L2Norm();
            if(d > 1)
            {
                return new double[0];
            }
            double psi = Math.Asin(d);
            
            return new double[2] { -phi+psi,- phi-psi+Math.PI };

        }

    }
    public enum JointType
    {
        Revolute,
        Prismatic,
        MobileOrientation,
        MobileTranslation
    }

    public class Robot
    {
        public Robot() { }
        public Matrix<double> H { get; set; }
        public Matrix<double> P { get; set; }
        public JointType[] Joint_type { get; set; }
        public double[] Joint_lower_limit { get; set; }
        public double[] Joint_upper_limit { get; set; }
        public double[] Joint_vel_limit { get; set; }
        public double[] Joint_acc_limit { get; set; }
        public Matrix<double>[] M { get; set; }
        public Matrix<double> R_tool { get; set; }
        public Vector<double> P_tool { get; set; }
        public string[] Joint_names { get; set; }
        public string Root_link_name { get; set; }
        public string Tip_link_name { get; set; }

        public Robot(Matrix<double> h, Matrix<double> p, JointType[] joint_type, double[] joint_lower_limit = default(double[]), double[] joint_upper_limit = default(double[]), double[] joint_vel_limit = default(double[]), double[] joint_acc_limit = default(double[]), Matrix<double>[] m = default(Matrix<double>[]),
                 Matrix<double> r_tool = default(Matrix<double>), Vector<double> p_tool = default(Vector<double>), string[] joint_names = default(string[]), string root_link_name = default(string), string tip_link_name = default(string))
        {
            IEnumerable<Vector<double>> splicer = h.EnumerateColumns();
            
            foreach (Vector<double> column in splicer)
            {
                if (!(column.L2Norm().AlmostEqual(1,8))) throw new ArgumentException(String.Format("Matrix H is not Acceptable"));
            }
            /*for (int i = 0; i < joint_type.Length; i++)
            {

                if (!(joint_types.Contains(joint_type[i]))) throw new ArgumentException(String.Format("Joint types contains incorrect values"));
            }*/
            if (h.RowCount != 3) throw new ArgumentException(String.Format("Matrix H is not Acceptable"));
            if (p.RowCount != 3) throw new ArgumentException(String.Format("Matrix P is not Acceptable"));
            if (h.ColumnCount + 1 != p.ColumnCount || h.ColumnCount != joint_type.Length) throw new ArgumentException(String.Format("Matrix Dimensions are not Acceptable"));
            if (joint_lower_limit != null && joint_upper_limit != null)
            {
                if (joint_lower_limit.Length != joint_type.Length || joint_type.Length != joint_upper_limit.Length) throw new ArgumentException(String.Format("Joint Limits not Acceptable"));
                Joint_upper_limit = joint_upper_limit;
                Joint_lower_limit = joint_lower_limit;
            }
            else
            {
                Joint_lower_limit = null;
                Joint_upper_limit = null;
            }
            if (joint_vel_limit != null)
            {
                if (joint_vel_limit.Length != joint_type.Length) throw new ArgumentException(String.Format("Joint Velocities not Acceptable"));
                Joint_vel_limit = joint_vel_limit;
            }
            else
            {
                Joint_vel_limit = null;
            }
            if (joint_acc_limit != null)
            {
                if (joint_acc_limit.Length != joint_type.Length) throw new ArgumentException(String.Format("Joint Accelerations not Acceptable"));
                Joint_acc_limit = joint_acc_limit;
            }
            else
            {
                Joint_acc_limit = null;
            }
            if (m != null)
            {
                if (m.Length != H.ColumnCount) throw new ArgumentException(String.Format("Inertia Matrices not Acceptable"));
                for (int i = 0; i < m.Length; i++)
                {
                    if (m[i].ColumnCount != 6 || m[i].RowCount != 6) throw new ArgumentException(String.Format("Inertia Matrices not Acceptable"));
                }
                M = m;
            }
            else
            {
                M = null;
            }
            if (r_tool != null && p_tool != null)
            {
                R_tool = r_tool;
                P_tool = p_tool;
            }
            else
            {
                R_tool = null;
                P_tool = null;
            }
            H = h;
            P = p;
            Joint_type = joint_type;
            if (joint_names != null)
            {
                if (joint_names.Length != joint_type.Length) throw new ArgumentException(String.Format("Joint Names not Acceptable"));
                Joint_names = joint_names;
            }
            else
            {
                Joint_names = null;
            }
            Root_link_name = root_link_name;
            Tip_link_name = tip_link_name;


        }

    }
    public class Transform
    {
        public Transform() { }
        public Matrix<double> R { get; set; }
        public Vector<double> P { get; set; }
        public string Parent_frame_id { get; set; }
        public string Child_frame_id { get; set; }

        public Transform(Matrix<double> r, Vector<double> p, string parent_frame_id = default(string), string child_frame_id = default(string))
        {
            if (r.ColumnCount != 3 || r.RowCount != 3) throw new ArgumentException(String.Format("Rotation Matrix for transform not Acceptable"));
            if(p.Count!=3) throw new ArgumentException(String.Format("Position Vector for transform not Acceptable"));
            R = r;
            P = p;
            Parent_frame_id = parent_frame_id;
            Child_frame_id = child_frame_id;
        }
        public static Transform operator *(Transform tran1, Transform tran2)
        {
            Matrix<double> R = tran1.R* tran2.R;
            Vector<double> P = tran1.P + tran1.R * tran2.P;
            Transform output= new Transform(R, P, tran1.Parent_frame_id, tran2.Child_frame_id);
            return output;
        }
        public static bool operator ==(Transform tran1, Transform tran2)
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; i < 3; i++)
                {
                    if (!(tran1.R[i, j].AlmostEqual(tran2.R[i, j], 8))) return false;
                }
            }
            for (int i = 0; i < 3; i++)
            {
                if (!(tran1.P[i].AlmostEqual(tran2.P[i]))) return false;
            }
            return true;

        }
        public static bool operator !=(Transform tran1, Transform tran2)
        {
            return (!(tran1 == tran2));
        }
        public Transform inv(Transform tran1)
        {
            Matrix<double> R = tran1.R.Transpose();
            Vector<double> P = -(tran1.R*tran1.P);
            Transform output = new Transform(R, P, tran1.Parent_frame_id, tran1.Child_frame_id);
            return output;
        }
    }
}
