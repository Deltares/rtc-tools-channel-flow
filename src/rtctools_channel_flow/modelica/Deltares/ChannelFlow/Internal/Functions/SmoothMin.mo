within Deltares.ChannelFlow.Internal.Functions;

function SmoothMin "Smooth Approximation of a Min() Function"
  input Real a;
  input Real b;
  // A small value to ensure smoothness
  input Real eps = Deltares.Constants.eps;
  output Real smooth_min;
algorithm
  smooth_min := -1.0 * Deltares.Functions.SmoothMax(a=-1.0 * a, b=-1.0 * b, eps=eps);
end SmoothMin;