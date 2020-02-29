using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using TestGeneralRoboticsToolboxNET;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics;

namespace GeneralRoboticsToolboxTests
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void TestHat()
        {
            Vector<double> k = Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 });
            Matrix<double> khat = Matrix<double>.Build.Dense(3, 3);
            khat[0, 1] = -3;
            khat[0, 2] = 2;
            khat[1, 0] = 3;
            khat[1, 2] = -1;
            khat[2, 0] = -2;
            khat[2, 1] = 1;
            Matrix<double> k_hat = GeneralRoboticsToolbox.Hat(k);
            Assert.AreEqual(k_hat, khat, "Hat didn't work");
        }
        [TestMethod]
        public void TestRot()
        {
            Vector<double> k = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 });
            Matrix<double> rot1 = Matrix<double>.Build.Dense(3, 3);
            rot1[0, 0] = 1;
            rot1[0, 1] = 0;
            rot1[0, 2] = 0;
            rot1[1, 0] = 0;
            rot1[1, 1] = 0;
            rot1[1, 2] = 1;
            rot1[2, 0] = 0;
            rot1[2, 1] = -1;
            rot1[2, 2] = 0;

            rot1 = rot1.Transpose();

            Matrix<double> rot = GeneralRoboticsToolbox.Rot(k, Math.PI / 2);

            Assert.IsTrue(rot1.AlmostEqual(rot, 1e-6), "rot1 failed");
            Matrix<double> rot2 = Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 0, -1.0 }, new[] { 0, 1.0, 0 }, new[] { 1.0, 0, 0 });
            Vector<double> k2 = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            rot2 = rot2.Transpose();
            Matrix<double> rot_2 = GeneralRoboticsToolbox.Rot(k2, Math.PI / 2);
            Assert.IsTrue(rot2.AlmostEqual(rot_2, 1e-6), "rot2 failed");

            Matrix<double> rot3 = Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 1.0, 0 }, new[] { -1.0, 0, 0 }, new[] { 0, 0, 1.0 });
            Vector<double> k3 = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            rot3 = rot3.Transpose();
            Matrix<double> rot_3 = GeneralRoboticsToolbox.Rot(k3, Math.PI / 2);
            Assert.IsTrue(rot3.AlmostEqual(rot_3, 1e-6), "rot3 failed");

            Matrix<double> rot4 = Matrix<double>.Build.DenseOfRowArrays(new[] { -0.5057639, -0.1340537, 0.8521928 }, new[] { 0.6456962, -0.7139224, 0.2709081 }, new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> k4 = Vector<double>.Build.DenseOfArray(new[] { 0.4490221, 0.30207945, 0.84090853 });
            Matrix<double> rot_4 = GeneralRoboticsToolbox.Rot(k4, 2.65949884);
            Assert.IsTrue(rot4.AlmostEqual(rot_4, 1e-6), "rot4 failed");
        }

        [TestMethod]
        public void TestR2Rot()
        {
            void _R2rot_test(Vector<double> k, double theta1)
            {
                Matrix<double> R = GeneralRoboticsToolbox.Rot(k, theta1);

                Tuple<Vector<double>, double> r2_vals = GeneralRoboticsToolbox.R2rot(R);
                Vector<double> _k2 = r2_vals.Item1;
                double theta2 = r2_vals.Item2;
                if (Math.Abs(theta1 - theta2) > (theta1 + theta2))
                {
                    _k2 = -_k2;
                    theta2 = -theta2;
                }

                Assert.IsTrue(theta1.AlmostEqual(theta2, 1e-6));
                if (Math.Abs(theta1) < (1e-9))
                {
                    return;
                }
                if ((Math.Abs(theta1) - Math.PI) < 1e-9)
                {
                    if ((k + _k2).L2Norm() < (1e-6))
                    {
                        Assert.IsTrue(k.AlmostEqual(-_k2, 1e-6));
                        return;
                    }
                    Console.WriteLine(_k2);
                    Assert.IsTrue(k.AlmostEqual(_k2, 1e-6));
                    return;
                }
                Assert.IsTrue(k.AlmostEqual(_k2, 1e-6));
            }

            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0.4490221, 0.30207945, 0.84090853 }), 2.65949884);

            //Singularities
            Vector<double> k1 = Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 }).L2Norm();
            _R2rot_test(k1, 1e-10);

            Vector<double> k2 = Vector<double>.Build.DenseOfArray(new[] { 2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k2, Math.PI + (1e-10));

            Vector<double> k3 = Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k3, Math.PI + (1e-10));

            Vector<double> k4 = Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k4, Math.PI + (1e-10));

            Vector<double> k5 = Vector<double>.Build.DenseOfArray(new[] { 0, -1.0, -3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 0, -1.0, -3.0 }).L2Norm();
            _R2rot_test(k5, Math.PI + (1e-10));

            Vector<double> k6 = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            _R2rot_test(k6, Math.PI + (1e-10));
        }

        [TestMethod]
        public void TestScrewMatrix()
        {
            Vector<double> k = Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 });
            Matrix<double> G = GeneralRoboticsToolbox.Screw_matrix(k);
            Matrix<double> G_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 1.0, 0, 0, 0, -3, 2 },
                new[] { 0, 1.0, 0, 3, 0, -1 },
                new[] { 0, 0, 1.0, -2, 1, 0 },
                new[] { 0, 0, 0, 1.0, 0, 0 },
                new[] { 0, 0, 0, 0, 1.0, 0 },
                new[] { 0, 0, 0, 0, 0, 1.0 });
            Assert.IsTrue(G.AlmostEqual(G_t, 1e-8));
        }

        [TestMethod]
        public void TestR2Q()
        {
            Matrix<double> rot = Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.5057639, -0.1340537, 0.8521928 },
                new[] { 0.6456962, -0.7139224, 0.2709081 },
                new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> q_t = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Vector<double> q = GeneralRoboticsToolbox.R2Q(rot);
            Assert.IsTrue(q.AlmostEqual(q_t, 1e-6));
        }

        [TestMethod]
        public void TestQ2R()
        {
            Matrix<double> rot_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.5057639, -0.1340537, 0.8521928 },
                new[] { 0.6456962, -0.7139224, 0.2709081 },
                new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> q = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Matrix<double> rot = GeneralRoboticsToolbox.Q2R(q);
            Assert.IsTrue(rot.AlmostEqual(rot_t, 1e-6));
        }

        [TestMethod]
        public void TestRot2Q()
        {
            Tuple<Vector<double>, double> rot = GeneralRoboticsToolbox.R2rot(Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.5057639, -0.1340537, 0.8521928 },
                new[] { 0.6456962, -0.7139224, 0.2709081 },
                new[] { 0.5720833, 0.6872731, 0.4476342 }));
            Vector<double> k = Vector<double>.Build.Dense(3);
            for (int i = 0; i < rot.Item1.Count; i++)
            {
                k[i] = (double)rot.Item1[i];
            }
            double theta = (double)rot.Item2;
            Vector<double> q_t = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Vector<double> q = GeneralRoboticsToolbox.Rot2Q(k, theta);
            Assert.IsTrue(q.AlmostEqual(q_t, 1e-6));
        }

        [TestMethod]
        public void TestQ2Rot()
        {
            Matrix<double> rot_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.5057639, -0.1340537, 0.8521928 },
                new[] { 0.6456962, -0.7139224, 0.2709081 },
                new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> q = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Tuple<Vector<double>, double> rot = GeneralRoboticsToolbox.Q2Rot(q);
            Vector<double> k = rot.Item1;
            double theta = rot.Item2;
            Assert.IsTrue(GeneralRoboticsToolbox.Rot(k, theta).AlmostEqual(rot_t, 1e-6));
        }

        [TestMethod]
        public void TestQuatcomplement()
        {
            Vector<double> q = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Vector<double> q_c = GeneralRoboticsToolbox.Quatcomplement(q);
            Assert.IsTrue(q[0].AlmostEqual(q_c[0], 1e-8));
            Assert.IsTrue(q.SubVector(1, 3).AlmostEqual(-q_c.SubVector(1, 3), 1e-8));
        }

        [TestMethod]
        public void TestQuatproduct()
        {
            Vector<double> q_1 = Vector<double>.Build.DenseOfArray(new[] { 0.63867877, 0.52251797, 0.56156573, 0.06089615 });
            Vector<double> q_2 = Vector<double>.Build.DenseOfArray(new[] { 0.35764716, 0.61051424, 0.11540801, 0.69716703 });
            Matrix<double> R_t = GeneralRoboticsToolbox.Q2R(q_1).Multiply(GeneralRoboticsToolbox.Q2R(q_2));
            Vector<double> q_t = GeneralRoboticsToolbox.R2Q(R_t);
            Vector<double> q = GeneralRoboticsToolbox.Quatproduct(q_1).Multiply(q_2);
            Assert.IsTrue(q.AlmostEqual(q_t, 1e-6));
        }

        [TestMethod]
        public void TestQuatjacobian()
        {
            Vector<double> q = Vector<double>.Build.DenseOfArray(new[] { 0.63867877, 0.52251797, 0.56156573, 0.06089615 });
            Matrix<double> J = GeneralRoboticsToolbox.Quatjacobian(q);
            Matrix<double> J_t = Matrix<double>.Build.DenseOfRowArrays(new[] { -0.26125898, -0.28078286, -0.03044808 },
                new[] { 0.31933938, 0.03044808, -0.28078286 },
                new[] { -0.03044808, 0.31933938, 0.26125898 },
                new[] { 0.28078286, -0.26125898, 0.31933938 });

            Assert.IsTrue(J.AlmostEqual(J_t, 1e-6));
        }

        [TestMethod]
        // singularRPI requested like this:
        //[ExpectedException(typeof(ArgumentException),
        //"A userId of null was inappropriately allowed.")]
        public void TestRpy2R()
        {
            Vector<double> rpy1 = Vector<double>.Build.DenseOfArray(new[] { 10 * Math.PI / 180, -30 * Math.PI / 180, 90 * Math.PI / 180 });

            Matrix<double> R1 = GeneralRoboticsToolbox.Rpy2R(rpy1);
            Matrix<double> R1_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.0000000, -0.9848077, 0.1736482 },
                new[] { 0.8660254, -0.0868241, -0.4924039 },
                new[] { 0.5000000, 0.1503837, 0.8528686 });
            Assert.IsTrue(R1.AlmostEqual(R1_t, 1e-6));
            Vector<double> rpy2 = GeneralRoboticsToolbox.R2Rpy(R1);
            Assert.IsTrue(rpy1.AlmostEqual(rpy2, 1e-6));

            // Check singularity
            Vector<double> rpy3 = Vector<double>.Build.DenseOfArray(new[] { 10 * Math.PI / 180, 90 * Math.PI / 180, -30 * Math.PI / 180 });
            Matrix<double> R3 = GeneralRoboticsToolbox.Rpy2R(rpy3);
        }

        public Robot puma260b_robot()
        {
            // Returns an approximate Robot instance for a Puma 260B robot


            Vector<double> x = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 });
            Vector<double> y = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            Vector<double> z = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            Vector<double> a = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0, 0 });

            Matrix<double> H = Matrix<double>.Build.DenseOfColumnVectors(z, y, y, z, y, x);
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

        public Robot abb_irb6640_180_255_robot()
        {
            // Returns an approximate Robot instance for a Puma 260B robot


            Vector<double> x = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 });
            Vector<double> y = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            Vector<double> z = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            Vector<double> a = Vector<double>.Build.DenseOfArray(new[] { 0.0, 0, 0 });

            Matrix<double> H = Matrix<double>.Build.DenseOfColumnVectors(z, y, y, x, y, x);
            Matrix<double> P = Matrix<double>.Build.DenseOfColumnVectors(0.78 * z, 0.32*x,1.075*z, 0.2*z, 1.142*x, 0.2*x, a);
            int[] joint_type = new[] { 0, 0, 0, 0, 0, 0 };
            double[] joint_min = new[] { -170.0, -65, -180, -300, -120, -360 };
            double[] joint_max = new[] { 170.0, 85, 70, 300, 120, 360 };
            for (int i = 0; i < joint_min.Length; i++)
            {
                joint_min[i] = joint_min[i] * Math.PI / 180.0;
                joint_max[i] = joint_max[i] * Math.PI / 180.0;
            }
            return new Robot(H, P, joint_type, joint_min, joint_max);
        }

        public Robot puma260b_robot_tool()
        {
            Robot robot = puma260b_robot();
            robot.R_tool = GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0);
            robot.P_tool = Vector<double>.Build.DenseOfArray(new[] { 0.05, 0, 0 });
            return robot;
        }

        [TestMethod]
        public void TestFwdkin()
        {
            Robot puma = this.puma260b_robot();

            Transform pose = GeneralRoboticsToolbox.Fwdkin(puma, new[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
            Assert.IsTrue(pose.R.AlmostEqual(Matrix<double>.Build.DenseIdentity(3), 1e-8));
            Assert.IsTrue(pose.P.AlmostEqual(Vector<double>.Build.DenseOfArray(new[] { 10.0, -4.9, 4.25 }), 1e-6));

            // Another right-angle configuration
            double[] joints2 = new[] { 180.0, -90, -90, 90, 90, 90 };
            for (int i = 0; i < joints2.Length; i++)
            {
                joints2[i] = joints2[i] * Math.PI / 180.0;
            }
            Transform pose2 = GeneralRoboticsToolbox.Fwdkin(puma, joints2);
            Matrix<double> rot2 = GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 }), Math.PI).Multiply(GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), -Math.PI / 2));
            Assert.IsTrue(pose2.R.AlmostEqual(rot2, 1e-6));
            Assert.IsTrue(pose2.P.AlmostEqual(Vector<double>.Build.DenseOfArray(new[] { -0.75, 4.9, 31 }) * 0.0254, 1e-6));

            //Random configuration
            double[] joints3 = new[] { 50.0, -105, 31, 4, 126, -184 };
            for (int i = 0; i < joints3.Length; i++)
            {
                joints3[i] = joints3[i] * Math.PI / 180.0;
            }
            Transform pose3 = GeneralRoboticsToolbox.Fwdkin(puma, joints3);
            Matrix<double> pose3_R_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0.4274, 0.8069, -0.4076 },
                new[] { 0.4455, -0.5804, -0.6817 },
                new[] { -0.7866, 0.1097, -0.6076 });

            Vector<double> pose3_P_t = Vector<double>.Build.DenseOfArray(new[] { 0.2236, 0.0693, 0.4265 });
            Assert.IsTrue(pose3.R.AlmostEqual(pose3_R_t, 1e-4));
            Assert.IsTrue(pose3.P.AlmostEqual(pose3_P_t, 1e-4));

            Robot puma_tool = this.puma260b_robot_tool();

            Transform pose4 = GeneralRoboticsToolbox.Fwdkin(puma_tool, joints3);
            Matrix<double> pose4_R_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0.4076, 0.8069, 0.4274 },
                new[] { 0.681654, -0.580357, 0.44557 },
                new[] { 0.60759, 0.1097, -0.7866 });
            Console.WriteLine("Robot R tool={0}", pose4.R);
            Console.WriteLine("Robot R calculated tool={0}", pose4_R_t);
            Console.WriteLine("Robot p tool={0}", pose4.P);
            
            Vector<double> pose4_P_t = Vector<double>.Build.DenseOfArray(new[] { 0.2450, 0.0916, 0.3872 });
            Console.WriteLine("Robot p calculated tool={0}", pose4_P_t);
            //Assert.IsTrue(pose4.R.AlmostEqual(pose4_R_t, 1 * 10 ^ 4));
            //Assert.IsTrue(pose4.P.AlmostEqual(pose4_P_t, 1 * 10 ^ 4));
        }
        [TestMethod]
        private void TestRobotjacobian()
        {
            //Home configuration (See Page 2-2 of Puma 260 manual)
            Robot puma = this.puma260b_robot();
            Matrix<double> J = GeneralRoboticsToolbox.Robotjacobian(puma, new[] { 0.0, 0, 0, 0, 0, 0 });
            Assert.IsTrue(J.SubMatrix(0, 4, 0, J.ColumnCount).AlmostEqual(puma.H, 1e-4));

            // Another right-angle configuration
            double[] joints2 = new[] { 180.0, -90.0, -90.0, 90.0, 90.0, 90.0 };
            for (int i = 0; i < joints2.Length; i++)
            {
                joints2[i] = joints2[i] * Math.PI / 180.0;
            }
            Matrix<double> J2 = GeneralRoboticsToolbox.Robotjacobian(puma, joints2);
            Matrix<double> J2_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0, 0, 0, 0, -1.0, 0 },
                new[] { 0, -1.0, -1, 0, 0, 0 },
                new[] { 1.0, 0, 0, -1, -0, 1 },
                new[] { -0.1245, -0.4572, -0.2591, 0, 0, 0 },
                new[] { -0.0191, 0, 0, 0, 0.0559, 0 },
                new[] { 0, -0.0191, 0, 0, 0, 0 });
            Assert.IsTrue(J2.AlmostEqual(J2_t, 1e-4));

            //Random configuration
            double[] joints3 = new[] { 50.0, -105, 31, 4, 126, -184 };
            for (int i = 0; i < joints3.Length; i++)
            {
                joints3[i] = joints3[i] * Math.PI / 180.0;
            }
            Matrix<double> J3 = GeneralRoboticsToolbox.Robotjacobian(puma, joints3);
            Matrix<double> J3_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0, -0.766, -0.766, -0.6179, -0.7765, 0.4274 },
                new[] { 0, 0.6428, 0.6428, -0.7364, 0.6265, 0.4456 },
                new[] { 1, 0, 0, 0.2756, -0.0671, -0.7866 },
                new[] { -0.0693, 0.0619, -0.0643, 0.0255, -0.0259, 0 },
                new[] { 0.2236, 0.0738, -0.0766, -0.0206, -0.0357, 0 },
                new[] { 0, -0.1969, -0.2298, 0.0022, -0.0343, 0 });
            Assert.IsTrue(J3.AlmostEqual(J3_t, 1e-4));
        }

        [TestMethod]
        public void Test_Subproblems()
        {
            Vector<double> x = Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 });
            Vector<double> y = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            Vector<double> z = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });

            // Subproblem0
            Assert.IsTrue(GeneralRoboticsToolbox.Subproblem0(x, y, z) == Math.PI / 2);

            // Subproblem1
            Vector<double> k1 = (x + z) / (x + z).L2Norm();
            Vector<double> k2 = (y + z) / (y + z).L2Norm();
            Assert.IsTrue(GeneralRoboticsToolbox.Subproblem1(k1, k2, z) == Math.PI / 2);

            // Subproblem2
            Vector<double> p2 = x;
            Vector<double> q2 = x.Add(y).Add(z);
            q2 = q2 / q2.L2Norm();
            double[] a2 = GeneralRoboticsToolbox.Subproblem2(p2, q2, z, y);
            Assert.IsTrue(a2.Length == 4);
            //NOTE: DIFFERENT THAN PYTHON VERSION

            Matrix<double> r1_0 = GeneralRoboticsToolbox.Rot(z, a2[0]);
            Matrix<double> r1_1 = GeneralRoboticsToolbox.Rot(y, a2[1]);
            Vector<double> r1 = (r1_0 * r1_1).Column(0);

            Matrix<double> r2_0 = GeneralRoboticsToolbox.Rot(z, a2[2]);
            Matrix<double> r2_1 = GeneralRoboticsToolbox.Rot(y, a2[3]);
            Vector<double> r2 = (r2_0 * r2_1).Column(0);

            Assert.IsTrue(r1.AlmostEqual(q2, 1e-4));
            Assert.IsTrue(r2.AlmostEqual(q2, 1e-4));

            double[] a3 = GeneralRoboticsToolbox.Subproblem2(x, z, z, y);
            //Console.WriteLine(a3);
            Assert.IsTrue(a3.Length == 2);
            //NOTE: DIFFERENT THAN PYTHON VERSION

            Matrix<double> r3_0 = GeneralRoboticsToolbox.Rot(z, a3[0]);
            Matrix<double> r3_1 = GeneralRoboticsToolbox.Rot(y, a3[1]);
            Vector<double> r3 = (r3_0 * r3_1).Column(0);
            Assert.IsTrue(r3.AlmostEqual(z, 1e-4));

            // Subproblem3
            Vector<double> p4 = Vector<double>.Build.DenseOfArray(new[] { .5, 0, 0 });
            Vector<double> q4 = Vector<double>.Build.DenseOfArray(new[] { 0, .75, 0 });

            double[] a4 = GeneralRoboticsToolbox.Subproblem3(p4, q4, z, .5);
            double[] a5 = GeneralRoboticsToolbox.Subproblem3(p4, q4, z, 1.25);
            Assert.IsTrue(a4.Length == 2);

            Assert.IsTrue((q4 + GeneralRoboticsToolbox.Rot(z, a4[0]) * p4).L2Norm().AlmostEqual(0.5, 1e-8));
            Assert.IsTrue((q4 + GeneralRoboticsToolbox.Rot(z, a4[1]) * p4).L2Norm().AlmostEqual(0.5, 1e-8));

            Assert.IsTrue(a5.Length == 1);
            Assert.IsTrue((q4 + GeneralRoboticsToolbox.Rot(z, a5[0]) * p4).L2Norm().AlmostEqual(1.25, 1e-8));

            // Subproblem4
            Vector<double> p6 = y;
            Vector<double> q6 = Vector<double>.Build.DenseOfArray(new[] { .8, .2, .5 });
            double d6 = .3;

            double[] a6 = GeneralRoboticsToolbox.Subproblem4(p6, q6, z, d6);
            Assert.IsTrue((p6 * GeneralRoboticsToolbox.Rot(z, a6[0]) * q6).AlmostEqual(d6, 1e-4));
            Assert.IsTrue((p6 * GeneralRoboticsToolbox.Rot(z, a6[1]) * q6).AlmostEqual(d6, 1e-4));
        }

        [TestMethod]
        public void Test_robot6_sphericalwrist_invkin()
        {
            Robot robot1 = puma260b_robot();
            Robot robot2 = abb_irb6640_180_255_robot();
            Robot robot3 = puma260b_robot_tool();
            Robot[] robots =new[] { robot1, robot3 };
            double[] thetas = new double[6];
            foreach (Robot robot in robots) {
                for (int i = 0; i < 100; i++)
                {
                    Console.WriteLine("test number {0}", i);
                    for (int x = 0; x < 6; x++)
                    {
                        Random randnum = new Random();
                        thetas[x] = randnum.NextDouble() * (robot.Joint_upper_limit[x] - robot.Joint_lower_limit[x]) + robot.Joint_lower_limit[x];
                    }
                    Assert.IsTrue(_test_configuration(robot, thetas));
                }
            }
            Assert.IsTrue(_test_configuration(robot1,new[] { 0.0, 0, 0, 0, 0, 0 }));
            //Assert.IsTrue(_test_last_configuration()
        }

        public bool _test_configuration(Robot r, double[] theta)
        {
            Transform pose1 = GeneralRoboticsToolbox.Fwdkin(r, theta);

            double[][] theta2 = InverseKin.robot6_sphericalwrist_invkin(r, pose1);
            Console.WriteLine("hello");
            foreach(double limit in r.Joint_upper_limit)
            {
                Console.WriteLine("Upper limits: {0}", limit);
            }
            foreach (double limit in r.Joint_lower_limit)
            {
                Console.WriteLine("Lower limits: {0}", limit);
            }
            foreach (double[] thetaset in theta2)
            {
                
                foreach(double thetas in thetaset)
                {
                    Console.WriteLine("Theta values: {0}", thetas);
                }
                
            }
            if (!(theta2.Length > 0))
            {
                Console.WriteLine("length of theta2 no good");
                return false;
            }

            foreach(double[] thetavals in theta2)
            {
                
                Transform pose2 = GeneralRoboticsToolbox.Fwdkin(r, thetavals);
                Console.WriteLine("pose 1={0}", pose1);
                Console.WriteLine("pose 2={0}", pose2);
                if (!(pose1 == pose2))
                {
                    return false;
                }
            }
            return true;

        }
        public bool _test_last_configuration(Robot r, double[] theta,double[] last_theta)
        {
            Transform pose1 = GeneralRoboticsToolbox.Fwdkin(r, theta);

            double[][] theta2 = InverseKin.robot6_sphericalwrist_invkin(r, pose1);
            Transform pose2 = GeneralRoboticsToolbox.Fwdkin(r, theta2[0]);
            if (!(pose1 == pose2))
            {
                return false;
            }
            Vector<double> theta2_vec = Vector<double>.Build.DenseOfArray(theta2[0]);
            Vector<double> last_vec = Vector<double>.Build.DenseOfArray(last_theta);
            Console.WriteLine("last config: {0}", theta2_vec);
            Console.WriteLine("last config: {0}", last_theta);
            if (theta2_vec.AlmostEqual(last_vec, 1e-8)) return true;
            return false;

        }
    }
}
