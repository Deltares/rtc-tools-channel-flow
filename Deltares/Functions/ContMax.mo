within Deltares.Functions;

function ContMax "Continuous Approximation of a Max() Function"
  input Real a;
  input Real b;
  output Real contmax;
algorithm
  contmax := sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2;
end ContMax;