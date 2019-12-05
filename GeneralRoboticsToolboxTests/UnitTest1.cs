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
            Assert.AreEqual(k_hat, khat,  "Hat didn't work");
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
            rot1[2, 2] =0;
            

            
            rot1 = rot1.Transpose();
            
            
            Matrix<double> rot = GeneralRoboticsToolbox.Rot(k,Math.PI/2);
            
            
            Assert.IsTrue(rot1.AlmostEqual(rot, 1*10^-8),"rot1 failed") ;
                
            Matrix<double> rot2= Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 0,-1.0 }, new[] { 0, 1.0,0 },new[] { 1.0, 0, 0 });
            Vector<double> k2 = Vector<double>.Build.DenseOfArray(new[] { 0, 1.0, 0 });
            rot2 = rot2.Transpose();
            Matrix<double> rot_2 = GeneralRoboticsToolbox.Rot(k2, Math.PI/2);
            Assert.IsTrue(rot2.AlmostEqual(rot_2, 1*10^-8), "rot2 failed");



            Matrix<double> rot3 = Matrix<double>.Build.DenseOfRowArrays(new[] { 0, 1.0, 0 }, new[] { -1.0, 0, 0 }, new[] { 0, 0, 1.0 });
            Vector<double> k3 = Vector<double>.Build.DenseOfArray(new[] { 0, 0, 1.0 });
            rot3 = rot3.Transpose();
            Matrix<double> rot_3 = GeneralRoboticsToolbox.Rot(k3, Math.PI / 2);

            Assert.IsTrue(rot3.AlmostEqual(rot_3, 1*10^-8), "rot3 failed");
            Matrix<double> rot4 = Matrix<double>.Build.DenseOfRowArrays(new[] { -0.5057639, -0.1340537, 0.8521928 }, new[] { 0.6456962, -0.7139224, 0.2709081 }, new[] { 0.5720833, 0.6872731, 0.4476342 });
            Vector<double> k4 = Vector<double>.Build.DenseOfArray(new[] { 0.4490221, 0.30207945, 0.84090853 });
            rot4 = rot4.Transpose();
            Matrix<double> rot_4 = GeneralRoboticsToolbox.Rot(k4, 2.65949884);
            Assert.IsTrue(rot4.AlmostEqual(rot_4, 1*10^-8), "rot4 failed");
        }

    

    }
}
