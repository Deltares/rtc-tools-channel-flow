within Deltares.ChannelFlow.Internal.Functions;

function SmoothAbs "Smooth Approximation of an Abs() Function"
  input Real a;
  // A small value to ensure smoothness
  input Real eps = Deltares.Constants.eps;
  output Real smooth_abs;
algorithm
  smooth_abs := sqrt(a ^ 2 + eps);
end SmoothAbs;