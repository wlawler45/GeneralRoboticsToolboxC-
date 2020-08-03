using GeneralRoboticsToolbox;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.Optimization;
using System;
using System.Collections.Generic;
using System.Text;

namespace GeneralRoboticsToolboxNET
{
    public class IterativeInverseKinSolver
    {
        public IterativeInverseKinSolver() { }

        

        public Robot Robot { get; set; }
        public IterativeInverseKinSolver(Robot robot)
        {
            Robot = robot;
        }

        static bool check_converge(double val)
        {
            return (Math.Abs(val) < 0.0001);
        }
        public  Tuple<double[], bool> Calculateinversekin(Transform desired_pose, double[] current_joint_angles, double epsilon = 0.001, int max_steps = 200, double step_size_alpha = 1.0)
        {

            int num_joints = Robot.Joint_type.Length;
            if (current_joint_angles.Length != num_joints)
            {
                Console.WriteLine("Current joint angles length is not valid");
            }

            Matrix<double> Kq = Matrix<Double>.Build.DenseIdentity(num_joints) * epsilon;
            int iteration_count = 0;
            bool converged = false;
            Vector<double> q_cur = Vector<double>.Build.DenseOfArray(current_joint_angles);
            while(iteration_count<max_steps && !(converged))
            {
                Transform pose = Functions.Fwdkin(Robot, q_cur.ToArray());
                Matrix<double> R_cur = pose.R;
                Vector<double> p_cur = pose.P;
                Matrix<double> J0T = Functions.Robotjacobian(Robot, q_cur.ToArray());
                Matrix<double> Tr = Matrix<double>.Build.Dense(6,6);
                Tr.SetSubMatrix(0, 0, R_cur.Transpose());
                Tr.SetSubMatrix(3, 3, R_cur.Transpose());
                J0T = Tr.Multiply(J0T);
                Matrix<double> ER = desired_pose.R.Transpose().Multiply(R_cur);
                Vector<double> EPfirststep = p_cur - desired_pose.P;
                Matrix<double> EPsecondstep = R_cur.Transpose();
                Vector<double> EP = EPsecondstep.Multiply(EPfirststep);
                Tuple<Vector<double>, double> returned = Functions.R2rot(ER);
                Vector<double> s = returned.Item1 * Math.Sin(returned.Item2 / 2);
                Matrix<double> deltastep1 = Kq + J0T.Transpose().Multiply(J0T);
                //Vector<double>[] array = [s, EP];
                Console.WriteLine(s.Count);
                Console.WriteLine(EP.Count);

                Matrix<double> deltastep2 = Matrix<double>.Build.DenseOfRowArrays(
                new[] { s[0], s[1], s[2], EP[0], EP[1], EP[2] });
                
                Matrix<double> deltastep3 = deltastep1.Inverse().Multiply(J0T.Transpose());
                Matrix<double> delta= step_size_alpha*(deltastep3.Multiply(deltastep2.Transpose()));
                Vector<double> deltavector = delta.Column(0);
                q_cur = q_cur - delta.Column(0);
                converged = deltastep2.ForAll(check_converge);
                iteration_count++;

            }
           

            return Tuple.Create(q_cur.ToArray(), converged);
        }

    }
}
