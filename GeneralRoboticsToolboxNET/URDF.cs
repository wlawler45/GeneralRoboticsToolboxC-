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
using System.Xml;
using System.Xml.Linq;
using System.Xml.XPath;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace TestGeneralRoboticsToolboxNET
{

    
    public struct Joint
    {
        public string Name;
        public Int32 Joint_type;
        public Transform Origin;
        public string Parent;
        public string Child;
        public Vector<Double> Axis;
        public double Effort_limit;
        public double Lower_limit;
        public double Upper_limit;
        public double Velocity_limit;

        /*public Joint(string joint_type,Transform origin,string parent,string child, Vector<double> axis,double effort_limit, double lower_limit, double upper_limit, double velocity_limit)
        {
            Joint_type=joint_type;
            Origin=origin;
            Parent=parent;
            Child=child;
            Axis=axis;
            Effort_limit=effort_limit;
            Lower_limit=lower_limit;
            Upper_limit=upper_limit;
            Velocity_limit=velocity_limit;
        }*/
    }
    public struct Link
    {
        public string Name;
        public Transform Origin;
        public double Mass;
        public double[] Inertia;
        public Matrix<double> InertiaMatrix;
        /*
        public Link(Transform origin,double mass, double[] inertia)
        {
            Origin = origin;
            Mass = mass;
            Inertia = inertia;
        }*/
    }
    class URDF
    {
        public XmlDocument _load_xml_file(string filename, string package = default(string))
        {
            XmlDocument doc = new XmlDocument();
            doc.PreserveWhitespace = true;
            try { doc.Load(filename); }
            catch(System.IO.FileNotFoundException)
            {
                throw new ArgumentException(String.Format("File name not found"));
            }
            return doc;
        }
        public List<string> find_possible_tip_link_names(Dictionary<string,Joint> joints)
        {
            HashSet<string> tip_link_candidate = new HashSet<string>();
            HashSet<string> tip_link_remove = new HashSet<string>();
            List<string> returns= new List<string>();
            foreach (KeyValuePair<string, Joint> i in joints)
            {
                tip_link_candidate.Add(i.Value.Child);
            }
            foreach (KeyValuePair<string, Joint> i in joints)
            {
                tip_link_candidate.Remove(i.Value.Parent);
            }
            foreach (string j in tip_link_candidate)
            {
                returns.Add(j);
            }
            return returns;

        }
        public List<string> find_possible_root_link_names(Dictionary<string,Joint> joints)
        {
            HashSet<string> root_link_candidate = new HashSet<string>();
            HashSet<string> root_link_remove = new HashSet<string>();
            List<string> returns = new List<string>();
            foreach (KeyValuePair<string, Joint> i in joints)
            {
                root_link_candidate.Add(i.Value.Parent);
            }
            foreach (KeyValuePair<string, Joint> i in joints)
            {
                root_link_candidate.Remove(i.Value.Child);
            }
            foreach (string j in root_link_candidate)
            {
                returns.Add(j);
            }
            return returns;

        }
        public Robot Xmldocument_to_robot(XmlDocument doc,string root_link=default(string),string tip_link=default(string))
        {
            Dictionary<string, Joint> joints = new Dictionary<string, Joint>();
            Dictionary<string, Link> links = new Dictionary<string, Link>();
            string current_link_name = null;
            string current_joint_name = null;
            string[] separator = { " " };
            Int32 count = 3;

            XmlReaderSettings settings = new XmlReaderSettings();
            using (XmlReader reader = XmlReader.Create("sawyer_base_urdf.xml", settings))
            {
                // reader.MoveToContent();

                XDocument xRoot = XDocument.Load(reader);
                // Console.WriteLine(xRoot.ToString());
                IEnumerable<XElement>linklist = xRoot.Descendants("link");
                foreach (XElement link in linklist)
                {
                    Link current_link = new Link();
                    current_link_name =(string)( link.Attribute("name"));
                    current_link.Name = current_link_name;
                    bool has_origin = false;
                    bool has_mass = false;
                    //Console.WriteLine(current_link_name);
                    if (link.HasElements)
                    {
                        if (link.XPathSelectElement("inertial") != null)
                        {

                            
                            XElement inertial = link.Element("inertial");
                            if (inertial.XPathSelectElement("origin") != null)
                            {
                                XElement origins = inertial.Element("origin");


                                string xyz = origins.Attribute("xyz").Value;



                                double[] values = Array.ConvertAll(xyz.Split(separator, count, StringSplitOptions.RemoveEmptyEntries), Double.Parse);
                                Vector<double> originxyz = Vector<double>.Build.DenseOfArray(values);
                                string rpy = origins.Attribute("rpy").Value;
                                double[] rpyvalues = Array.ConvertAll(xyz.Split(separator, count, StringSplitOptions.RemoveEmptyEntries), Double.Parse);
                                Vector<double> originrpy = Vector<double>.Build.DenseOfArray(rpyvalues);
                                Matrix<double> R = GeneralRoboticsToolbox.Rpy2R(originrpy);
                                Transform origin = new Transform(R, originxyz);
                                current_link.Origin = origin;
                                has_origin = true;
                            }
                            if (inertial.XPathSelectElement("mass") != null)
                            {
                                XElement mass = inertial.Element("mass");
                                string mass_val = mass.Attribute("value").Value;
                                current_link.Mass = Convert.ToDouble(mass_val);
                                has_mass = true;
                            }
                            if (inertial.XPathSelectElement("inertia") != null)
                            {
                                if (!(has_mass&&has_origin))
                                {
                                    throw new ArgumentException(String.Format("Inertia matrix requires value for mass and origin to generate for link {0}", current_link_name));
                                }
                                XElement inertially = inertial.Element("inertia");
                                double[] inertias = new double[6];
                                Matrix<double> Inertia= Matrix<double>.Build.Dense(6,6);
                                Inertia[0,0] = Convert.ToDouble(inertially.Attribute("ixx").Value);
                                Inertia[0,1] =- Convert.ToDouble(inertially.Attribute("ixy").Value);
                                Inertia[0,2] =- Convert.ToDouble(inertially.Attribute("ixz").Value);
                                Inertia[1,0] =- Convert.ToDouble(inertially.Attribute("ixy").Value);
                                Inertia[1,1] = Convert.ToDouble(inertially.Attribute("iyy").Value);
                                Inertia[1,2] = -Convert.ToDouble(inertially.Attribute("iyz").Value);
                                Inertia[2, 0] = -Convert.ToDouble(inertially.Attribute("ixz").Value);
                                Inertia[2, 1] = -Convert.ToDouble(inertially.Attribute("iyz").Value);
                                Inertia[2, 2] = Convert.ToDouble(inertially.Attribute("izz").Value);
                                //Matrix<double> cross = GeneralRoboticsToolbox.Cross(current_link.Origin.P, current_link.Origin.P);
                                
                                current_link.InertiaMatrix = Inertia;
                            }

                        }
                        
                        
                    }
                    links.Add(current_link_name, current_link);

                }
                IEnumerable<XElement> jointlist = xRoot.Descendants("joint");
                foreach (XElement joint in jointlist)
                {
                    Joint current_joint = new Joint();
                    current_joint_name = (string)(joint.Attribute("name"));
                    current_joint.Name = current_joint_name;
                    //Console.WriteLine(current_joint_name);
                    string type = (string)(joint.Attribute("type"));
                    switch (type)
                    {
                        case "revolute":
                            current_joint.Joint_type = 0;
                            break;
                        case "fixed":
                            current_joint.Joint_type = 4;
                            break;
                        default:
                            throw new ArgumentException(String.Format("Joint type for joint {0} not specified correctly", current_joint_name));

                    }
                    if (joint.HasElements)
                    {
                        if (joint.XPathSelectElement("origin") != null)
                        {
                            XElement origins = joint.Element("origin");


                            string xyz = origins.Attribute("xyz").Value;



                            double[] values = Array.ConvertAll(xyz.Split(separator, count, StringSplitOptions.RemoveEmptyEntries), Double.Parse);
                            Vector<double> originxyz = Vector<double>.Build.DenseOfArray(values);
                            string rpy = origins.Attribute("rpy").Value;
                            double[] rpyvalues = Array.ConvertAll(xyz.Split(separator, count, StringSplitOptions.RemoveEmptyEntries), Double.Parse);
                            Vector<double> originrpy = Vector<double>.Build.DenseOfArray(rpyvalues);
                            Matrix<double> R = GeneralRoboticsToolbox.Rpy2R(originrpy);
                            Transform origin = new Transform(R, originxyz);
                            current_joint.Origin = origin;

                        }
                        if (joint.XPathSelectElement("parent") != null)
                        {
                            XElement parent = joint.Element("parent");

                            current_joint.Parent = parent.Attribute("link").Value;

                        }
                        if (joint.XPathSelectElement("child") != null)
                        {
                            XElement child = joint.Element("child");

                            current_joint.Child = child.Attribute("link").Value;

                        }
                        if (joint.XPathSelectElement("limit") != null)
                        {
                            XElement limit = joint.Element("limit");
                            //Console.WriteLine("found limits!!!");
                            string effort = (string)limit.Attribute("effort");
                            if (effort != null) current_joint.Effort_limit = Convert.ToDouble(effort);

                            string lower = (string)limit.Attribute("lower");
                            if (lower != null) current_joint.Lower_limit = Convert.ToDouble(lower);
                            string upper = (string)limit.Attribute("upper");
                            if (upper != null) current_joint.Upper_limit = Convert.ToDouble(upper);
                            string velocity = (string)limit.Attribute("velocity");
                            if (velocity != null) current_joint.Velocity_limit = Convert.ToDouble(velocity);


                        }
                        if (joint.XPathSelectElement("axis") != null)
                        {
                            XElement axis = joint.Element("axis");
                            string xyz = axis.Attribute("xyz").Value;
                            double[] values = Array.ConvertAll(xyz.Split(separator, count, StringSplitOptions.RemoveEmptyEntries), Double.Parse);
                            Vector<double> axisxyz = Vector<double>.Build.DenseOfArray(values);
                            current_joint.Axis = axisxyz;
                        }
                    }
                    joints.Add(current_joint_name, current_joint);
                }
                

                string current_linked = root_link;
                bool root_not_specified = root_link == default(string);
                
                List<string> fails = new List<string>();
                Joint last_fail=new Joint();
                List<Joint> robot_joints = new List<Joint>();
                List<Joint> robot_joints_temp = new List<Joint>();
                if (tip_link != default(string))
                {
                    Console.WriteLine("Tip link found");
                    while (current_linked != tip_link)
                    {
                        foreach (KeyValuePair<string, Joint> linker in joints)
                        {
                            if ((linker.Value.Parent == current_linked || root_not_specified) && !fails.Contains(linker.Value.Name))
                            {
                                root_not_specified = false;
                                current_linked = linker.Value.Child;
                                
                                Console.WriteLine("name={0}", linker.Value.Name);
                                Console.WriteLine("searching for={0}", current_linked);
                                robot_joints.Add(linker.Value);
                                if (current_linked == tip_link)
                                {
                                    Console.WriteLine("finishing");
                                    break;
                                }


                                last_fail = linker.Value;



                                foreach (Joint joit in robot_joints)
                                {
                                    Console.WriteLine("current robot joints: {0}", joit.Name);
                                }




                            }
                        }
                        root_not_specified = root_link == default(string);


                        Console.WriteLine("adding to fails {0}", last_fail.Name);
                        fails.Add(last_fail.Name);
                        foreach (string j in fails)
                        {
                            Console.WriteLine("fails includes: {0}", j);
                        }

                        if (current_linked == tip_link)
                        {
                            Console.WriteLine("ending while");
                            break;
                        }
                        Console.WriteLine("escaped outer loop");

                        robot_joints.Clear();
                        current_linked = root_link;
                    }
                }
                else
                {
                    
                    List<string> possible_tips = find_possible_tip_link_names(joints);
                    foreach (string j in possible_tips)
                    {
                        Console.WriteLine("Possible tips: {0}", j);
                    }
                    
                       List<string> possible_roots = find_possible_root_link_names(joints);
                        foreach (string j in possible_roots)
                        {
                            Console.WriteLine("Possible roots: {0}", j);
                        }

                        foreach (string i in possible_roots)
                        {
                            current_linked = i;
                            bool any_answers = false;
                            string repeat = default(string);
                            while (!possible_tips.Contains(current_linked))
                            {
                                any_answers = false;
                                foreach (KeyValuePair<string, Joint> linker in joints)
                                {
                                    if ((linker.Value.Parent == current_linked) && !fails.Contains(linker.Value.Name))
                                    {
                                        any_answers = true;
                                        if (linker.Value.Joint_type == 4)
                                        {
                                            //Console.WriteLine("fixed joint in chain found: {0}", linker.Value.Name);
                                            last_fail = linker.Value;
                                            // fails.Add(last_fail.Name);

                                        }
                                        current_linked = linker.Value.Child;

                                        //Console.WriteLine("name={0}", linker.Value.Name);
                                        //Console.WriteLine("searching for={0}", current_linked);
                                        robot_joints.Add(linker.Value);
                                        if (possible_tips.Contains(current_linked))
                                        {
                                            //Console.WriteLine("finishing");
                                            break;
                                        }


                                        last_fail = linker.Value;



                                        




                                    }
                                }

                                if (!any_answers)
                                {
                                    //Console.WriteLine("no answers found for this loop");
                                    if (last_fail.Name == repeat)
                                    {
                                        goto Exit;
                                        
                                    }
                                    repeat = last_fail.Name;

                                }

                                //Console.WriteLine("adding to fails {0}", last_fail.Name);
                                fails.Add(last_fail.Name);
                                //foreach (string j in fails)
                                //{
                                    //Console.WriteLine("fails includes: {0}", j);
//                                }

                                if (possible_tips.Contains(current_linked))
                                {
                                    
                                    List<int> indices = new List<int>();
                                    
                                    for(int q=0;q<robot_joints.Count;q++)
                                    {
                                        if (robot_joints[q].Joint_type == 4)
                                        {
                                            indices.Add(q);
                                        }
                                    }
                                    if (indices.Count == 0)
                                    {
                                        if (robot_joints_temp.Count < robot_joints.Count)
                                        {
                                            
                                            robot_joints_temp = robot_joints;
                                        }
                                    }
                                    
                                    if (indices.Count ==1 && robot_joints.Count>2)
                                    {
                                        if(indices[0]==0 || indices[0] == robot_joints.Count - 1)
                                        {
                                            robot_joints.RemoveAt(indices[0]);
                                            
                                            if (robot_joints_temp.Count < robot_joints.Count)
                                            {
                                                
                                                robot_joints_temp = robot_joints;
                                            }
                                        }
                                    }
                                    if (indices.Count==2 && robot_joints.Count>4)
                                    {
                                        if(indices[0]==0 && indices[1] == (robot_joints.Count - 1))
                                        {
                                            robot_joints.RemoveAt(indices[0]);
                                            


                                            robot_joints.RemoveAt(indices[1]-1);
                                            
                                            foreach (Joint joit in robot_joints_temp)
                                            {
                                                
                                            }
                                            if (robot_joints_temp.Count < robot_joints.Count)
                                            {
                                                
                                                robot_joints_temp =new List<Joint> (robot_joints);
                                            }
                                        }
                                        
                                        
                                    }
                                }
                                
                                robot_joints.Clear();
                                
                                current_linked = i;
                            }
                        }
                    
                Exit:
                    robot_joints = robot_joints_temp;
                    
                    
                }

                foreach (Joint joit in robot_joints)
                {
                    Console.WriteLine("final robot joints: {0}", joit.Name);
                }
                Matrix<double> h = Matrix<double>.Build.Dense(3, robot_joints.Count);
                Matrix<double> p = Matrix<double>.Build.Dense(3, robot_joints.Count+1);
                Int32[] joint_types = new int[robot_joints.Count];
                double[] joint_lower_limits = new double[robot_joints.Count];
                double[] joint_upper_limits = new double[robot_joints.Count];
                double[] joint_vel_limits = new double[robot_joints.Count];
                double[] joint_acc_limits = new double[robot_joints.Count];
                string[] joint_names = new string[robot_joints.Count];
                
                for(int o = 0; o < robot_joints.Count; o++)
                {
                    h.SetColumn(o, robot_joints[o].Axis);
                    p.SetColumn(o, robot_joints[o].Origin.P);
                    joint_types[o] = robot_joints[o].Joint_type;
                    joint_lower_limits[o] = robot_joints[o].Lower_limit;
                    joint_upper_limits[o] = robot_joints[o].Upper_limit;
                    joint_vel_limits[o] = robot_joints[o].Velocity_limit;
                    joint_acc_limits[o] = robot_joints[o].Effort_limit;
                    joint_names[o] = robot_joints[o].Name;
                }
                
                string root_link_name;
                if (root_not_specified) root_link_name = robot_joints[0].Parent; else root_link_name = root_link;
                
                string tip_link_name;
                if (tip_link != default(string)) tip_link_name = robot_joints[-1].Child; else tip_link_name = tip_link;
                Robot robot = new Robot(h,p,joint_types,joint_lower_limits,joint_upper_limits,joint_vel_limits,joint_acc_limits);
                robot.Joint_names = joint_names;
                robot.Root_link_name = root_link_name;
                robot.Tip_link_name = tip_link_name;
                return robot;
            }


                
        }

    }
}
