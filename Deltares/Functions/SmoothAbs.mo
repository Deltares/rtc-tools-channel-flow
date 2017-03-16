within Deltares.Functions;

function SmoothAbs "Smooth Approximation of an Abs() Function"
  input Real a;
  output Real smoothabs;
algorithm
  smoothabs := sqrt(a^2 + Modelica.Constants.eps);
end SmoothAbs;