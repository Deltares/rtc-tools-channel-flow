within Deltares.ChannelFlow.Internal.Functions;

function SmoothMin "Smooth Approximation of a Min() Function"
  input Real a;
  input Real b;
  output Real smoothmin;
algorithm
  smoothmin := -1.0 * Deltares.Functions.SmoothMax( -1.0 * a, -1.0 * b);
end SmoothMin;