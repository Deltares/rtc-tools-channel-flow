within Deltares.Functions;

function ContMin "Continuous Approximation of a Min() Function"
  input Real a;
  input Real b;
  output Real contmin;
algorithm
  contmin := -1.0 * Deltares.Functions.ContMax( -1.0 * a, -1.0 * b);
end ContMin;