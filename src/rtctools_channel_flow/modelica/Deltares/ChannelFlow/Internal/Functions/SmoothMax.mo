within Deltares.ChannelFlow.Internal.Functions;

function SmoothMax "Smooth Approximation of a Max() Function"
  input Real a;
  input Real b;
  // A small value to ensure smoothness
  input Real eps = Deltares.Constants.eps;
  output Real smooth_max;
algorithm
  smooth_max := sqrt((a - b) ^ 2 + eps) / 2 + (a + b) / 2;
end SmoothMax;