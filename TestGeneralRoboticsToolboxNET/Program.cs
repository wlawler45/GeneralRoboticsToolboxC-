using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace TestGeneralRoboticsToolboxNET
{
    class Program
    {
        public static Robot puma260b_robot()
        {
            // Returns an approximate Robot instance for a Puma 260B robot


            Vector<double> x = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 });
            Vector<double> y = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            Vector<double> z = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            Vector<double> a = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0, 0 });

            Matrix<double> H = Matrix<double>.Build.DenseOfColumnVectors(x, y, y, z, y, x);
            Matrix<double> P = 0.0254 * Matrix<double>.Build.DenseOfColumnVectors(13 * z, a, (-4.9 * y + 7.8 * x - 0.75 * z), -8.0 * z, a, a, 2.2 * x);
            int[] joint_type = new[] { 0, 0, 0, 0, 0, 0 };
            double[] joint_min = new[] { -5.0, -256, -214, -384, -32, -267 };
            double[] joint_max = new[] { 313.0, 76, 34, 194, 212, 267 };
            for (int i = 0; i < joint_min.Length; i++)
            {
                joint_min[i] = joint_min[i] * Math.PI / 180.0;
                joint_max[i] = joint_max[i] * Math.PI / 180.0;
            }
            return new Robot(H, P, joint_type, joint_min, joint_max);
        }

        public static Robot puma260b_robot_tool()
        {
            Robot robot = puma260b_robot();
            robot.R_tool = GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0);
            robot.P_tool = Vector<double>.Build.DenseOfArray(new[] { 0.05, 0, 0 });
            return robot;
        }

        static void Main(string[] args)
        {
            var test = new GeneralRoboticsToolbox();

            double[] joints3 = new[] { 50.0, -105, 31, 4, 126, -184 };
            for (int i = 0; i < joints3.Length; i++)
            {
                joints3[i] = joints3[i] * Math.PI / 180.0;
            }

            Robot puma_tool = puma260b_robot_tool();

            Transform pose4 = GeneralRoboticsToolbox.Fwdkin(puma_tool, joints3);
            Matrix<double> pose4_R_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0.4076, 0.8069, 0.4274 },
                new[] { 0.6817, -0.5804, 0.4455 },
                new[] { 0.6076, 0.1097, -0.7866 });
            Vector<double> pose4_P_t = Vector<double>.Build.DenseOfArray(new[] { 0.2450, 0.0916, 0.3872 });
            //Console.WriteLine(GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0));
            Console.WriteLine(pose4.R);
            Console.ReadLine();
        }
    }
}


