within Deltares.Functions;

function SmoothMax "Smooth Approximation of a Max() Function"
  input Real a;
  input Real b;
  output Real smoothmax;
algorithm
  smoothmax := sqrt(((a) - (b)) ^ 2 + Modelica.Constants.eps) / 2 + ((a) + (b)) / 2;
end SmoothMax;