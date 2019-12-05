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
        static void Main(string[] args)
        {
            var test = new GeneralRoboticsToolbox();
            var urdf_test = new URDF();
            var doc=urdf_test._load_xml_file("sawyer_base_urdf.xml");
            urdf_test.Xmldocument_to_robot(doc, root_link: "right_arm_base_link", tip_link: "right_hand");
            Vector<double> k = Vector<double>.Build.DenseOfArray(new[] { 1.0, 2.0, 3.0}); 
            Matrix<double> k_hat = GeneralRoboticsToolbox.Hat(k);

            Console.WriteLine("{0}", k_hat);
            Console.ReadLine();
        }
    }
}


