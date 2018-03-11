within Deltares.ChannelFlow.Internal.Functions;

function SmoothSwitch
  input Real x;
  output Real y;
protected
  Real k = 50;
algorithm
  y := if x < -1 then 0 elseif x > 1 then 1 else 0 + (1 - 0) / (1 + exp(-k * x));
end SmoothSwitch;
