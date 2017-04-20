within Deltares.Functions;

function SmoothSwitch
  input Real a;
  output Real smoothswitch;
algorithm
  smoothswitch := if a < -1 then 0 elseif a > 1 then 1 else 0 + (1 - 0) / (1 + exp (Deltares.Functions.SmoothMin(-500 * a , 100 )));
end SmoothSwitch;
