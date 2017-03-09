within Deltares.Functions;

function ContAbs "Continuous Approximation of an Abs() Function"
  input Real a;
  output Real contabs;
algorithm
  contabs := sqrt(a^2 + Modelica.Constants.eps);
end ContAbs;