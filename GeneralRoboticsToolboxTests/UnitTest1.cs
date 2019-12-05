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
            rot1[2, 1] = 1;
            rot1[2, 2] = 0;



            rot1 = rot1.Transpose();


            Matrix<double> rot = GeneralRoboticsToolbox.Rot(k, Math.PI / 2);


            Assert.IsTrue(rot1.AlmostEqual(rot, 1 * 10 ^ -8), "rot1 failed");

            Matrix<double> rot2 = Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 0, -1.0 }, new[] { 0, 1.0, 0 }, new[] { 1.0, 0, 0 });
            Vector<double> k2 = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            rot2 = rot2.Transpose();
            Matrix<double> rot_2 = GeneralRoboticsToolbox.Rot(k2, Math.PI / 2);
            Assert.IsTrue(rot2.AlmostEqual(rot_2, 1 * 10 ^ -8), "rot2 failed");



            Matrix<double> rot3 = Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 1.0, 0 }, new[] { -1.0, 0, 0 }, new[] { 0, 0, 1.0 });
            Vector<double> k3 = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            rot3 = rot3.Transpose();
            Matrix<double> rot_3 = GeneralRoboticsToolbox.Rot(k3, Math.PI / 2);

            Assert.IsTrue(rot3.AlmostEqual(rot_3, 1 * 10 ^ -8), "rot3 failed");
            Matrix<double> rot4 = Matrix<double>.Build.DenseOfRowArrays(new[] { -0.5057639, -0.1340537, 0.8521928 }, new[] { 0.6456962, -0.7139224, 0.2709081 }, new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> k4 = Vector<double>.Build.DenseOfArray(new[] { 0.4490221, 0.30207945, 0.84090853 });
            rot4 = rot4.Transpose();
            Matrix<double> rot_4 = GeneralRoboticsToolbox.Rot(k4, 2.65949884);
            Assert.IsTrue(rot4.AlmostEqual(rot_4, 1 * 10 ^ -8), "rot4 failed");
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
                Assert.IsTrue(theta1.AlmostEqual(theta2, 1 * 10 ^ -6));
                if (Math.Abs(theta1) < (1 * 10 ^ -9))
                {
                    return;
                }
                if ((Math.Abs(theta1) - Math.PI) < (1 * 10 ^ -9))
                {
                    if ((k + _k2).L2Norm() < (1 * 10 ^ -9))
                    {
                        Assert.IsTrue(k.AlmostEqual(-_k2, 1 * 10 ^ -6));
                        return;
                    }
                    Assert.IsTrue(k.AlmostEqual(_k2, 1 * 10 ^ -6));
                    return;
                }
                Assert.IsTrue(k.AlmostEqual(_k2, 1 * 10 ^ -6));
            }

            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 1.0, 0, 0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 }), Math.PI / 2.0);
            _R2rot_test(Vector<double>.Build.DenseOfArray(new[] { 0.4490221, 0.30207945, 0.84090853 }), 2.65949884);

            //Singularities
            Vector<double> k1 = Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0 }).L2Norm();
            _R2rot_test(k1, 1 * 10 ^ -10);

            Vector<double> k2 = Vector<double>.Build.DenseOfArray(new[] { 2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k2, Math.PI + (1 * 10 ^ -10));

            Vector<double> k3 = Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k3, Math.PI + (1 * 10 ^ -10));

            Vector<double> k4 = Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }) / Vector<double>.Build.DenseOfArray(new[] { -2.0, -1.0, 3.0 }).L2Norm();
            _R2rot_test(k4, Math.PI + (1 * 10 ^ -10));

            Vector<double> k5 = Vector<double>.Build.DenseOfArray(new[] { 0, -1.0, -3.0 }) / Vector<double>.Build.DenseOfArray(new[] { 0, -1.0, -3.0 }).L2Norm();
            _R2rot_test(k5, Math.PI + (1 * 10 ^ -10));

            Vector<double> k6 = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            _R2rot_test(k6, Math.PI + (1 * 10 ^ -10));
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
            Assert.IsTrue(G.AlmostEqual(G_t, 1 * 10 ^ -8));
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
            Assert.IsTrue(q.AlmostEqual(q_t, 1 * 10 ^ -6));
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
            Assert.IsTrue(rot.AlmostEqual(rot_t, 1 * 10 ^ -6));
        }

        [TestMethod]
        public void TestRot2Q()
        {
            Tuple<Vector<double>, double> rot = GeneralRoboticsToolbox.R2rot(Matrix<double>.Build.DenseOfRowArrays(
                new[] { -0.5057639, -0.1340537, 0.8521928 },
                new[] { 0.6456962, -0.7139224, 0.2709081 },
                new[] { 0.5720833, 0.6872731, 0.4476342 }));
            float[] k = new float[3];
            for (int i = 0; i < rot.Item1.Count; i++)
            {
                k[i] = (float)rot.Item1[i];
            }
            float theta = (float)rot.Item2;
            Vector<float> q_t = Vector<float>.Build.DenseOfArray(new[] { (float)0.2387194, (float)0.4360402, (float)0.2933459, (float)0.8165967 });
            Vector<float> q = GeneralRoboticsToolbox.Rot2Q(k, theta);
            Assert.IsTrue(q.AlmostEqual(q_t, 1 * 10 ^ -6));
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
            Assert.IsTrue(GeneralRoboticsToolbox.Rot(k, theta).AlmostEqual(rot_t, 1 * 10 ^ -6));
        }

        [TestMethod]
        public void TestQuatcomplement()
        {
            Vector<double> q = Vector<double>.Build.DenseOfArray(new[] { 0.2387194, 0.4360402, 0.2933459, 0.8165967 });
            Vector<double> q_c = GeneralRoboticsToolbox.Quatcomplement(q);
            Assert.IsTrue(q[0].AlmostEqual(q_c[0], 1 * 10 ^ -8));
            Assert.IsTrue(q.SubVector(1, 3).AlmostEqual(-q_c.SubVector(1, 3), 1 * 10 ^ -8));
        }

        [TestMethod]
        public void TestQuatproduct()
        {
            Vector<double> q_1 = Vector<double>.Build.DenseOfArray(new[] { 0.63867877, 0.52251797, 0.56156573, 0.06089615 });
            Vector<double> q_2 = Vector<double>.Build.DenseOfArray(new[] { 0.35764716, 0.61051424, 0.11540801, 0.69716703 });
            Matrix<double> R_t = GeneralRoboticsToolbox.Q2R(q_1).Multiply(GeneralRoboticsToolbox.Q2R(q_2));
            Vector<double> q_t = GeneralRoboticsToolbox.R2Q(R_t);
            Vector<double> q = GeneralRoboticsToolbox.Quatproduct(q_1).Multiply(q_2);
            Assert.IsTrue(q.AlmostEqual(q_t, 1 * 10 ^ -6));
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

            Assert.IsTrue(J.AlmostEqual(J_t, 1 * 10 ^ -6));
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
            Assert.IsTrue(R1.AlmostEqual(R1_t, 1 * 10 ^ -6));
            Vector<double> rpy2 = GeneralRoboticsToolbox.R2Rpy(R1);
            Assert.IsTrue(rpy1.AlmostEqual(rpy2, 1 * 10 ^ -6));

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

        public Robot puma260b_robot_tool()
        {
            Robot robot = this.puma260b_robot();
            robot.R_tool = GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), Math.PI / 2.0);
            robot.P_tool = Vector<double>.Build.DenseOfArray(new[] { 0.05, 0, 0 });
            return robot;
        }

        [TestMethod]
        public void TestFwdkin()
        {
            Robot puma = this.puma260b_robot();

            Transform pose = GeneralRoboticsToolbox.Fwdkin(puma, new[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
            Assert.IsTrue(pose.R.AlmostEqual(Matrix<double>.Build.DenseIdentity(3), 1 * 10 ^ -8));
            Assert.IsTrue(pose.P.AlmostEqual(Vector<double>.Build.DenseOfArray(new[] { 10.0, -4.9, 4.25 }), 1 * 10 ^ -6));

            // Another right-angle configuration
            double[] joints2 = new[] { 180.0, -90, -90, 90, 90, 90 };
            for (int i = 0; i < joints2.Length; i++)
            {
                joints2[i] = joints2[i] * Math.PI / 180.0;
            }
            Transform pose2 = GeneralRoboticsToolbox.Fwdkin(puma, joints2);
            Matrix<double> rot2 = GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 }), Math.PI).Multiply(GeneralRoboticsToolbox.Rot(Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 }), -Math.PI / 2));
            Assert.IsTrue(pose2.R.AlmostEqual(rot2, 1 * 10 ^ -6));
            Assert.IsTrue(pose2.P.AlmostEqual(Vector<double>.Build.DenseOfArray(new[] { -0.75, 4.9, 31 }) * 0.0254, 1 * 10 ^ -6));

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
            Assert.IsTrue(pose3.R.AlmostEqual(pose3_R_t, 1 * 10 ^ -4));
            Assert.IsTrue(pose3.P.AlmostEqual(pose3_P_t, 1 * 10 ^ -4));

            Robot puma_tool = this.puma260b_robot_tool();

            Transform pose4 = GeneralRoboticsToolbox.Fwdkin(puma_tool, joints3);
            Matrix<double> pose4_R_t = Matrix<double>.Build.DenseOfRowArrays(
                new[] { 0.4076, 0.8069, 0.4274 },
                new[] { 0.6817, -0.5804, 0.4455 },
                new[] { 0.6076, 0.1097, -0.7866 });

            Vector<double> pose4_P_t = Vector<double>.Build.DenseOfArray(new[] { 0.2450, 0.0916, 0.3872 });
            Assert.IsTrue(pose4.R.AlmostEqual(pose4_R_t, 1 * 10 ^ 4));
            Assert.IsTrue(pose4.P.AlmostEqual(pose4_P_t, 1 * 10 ^ 4));
        }
    }
}
