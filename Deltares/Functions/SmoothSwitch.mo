within Deltares.Functions;

function SmoothSwitch
  input Real a;
  output Real smoothswitch;
algorithm
  smoothswitch := 0 + (1 - 0) / (1 + exp(-10 * a));
end SmoothSwitch;
