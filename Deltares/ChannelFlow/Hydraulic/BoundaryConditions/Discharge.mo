within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model Discharge "Defines a discharge"
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  function smoothswitch = Deltares.Functions.SmoothSwitch;
  input Modelica.SIunits.VolumeFlowRate Q;
  input Modelica.SIunits.MassFlowRate M[HQ.medium.n_substances];
  parameter Boolean upwind = true; // If true and there is outlfow from the system (into the discharge boudnary) then the concentration of the connected element is used.
equation
  HQ.Q + Q = 0;
  if upwind then
    HQ.M = -smoothswitch(Q) * M + (smoothswitch(Q) - 1) * Q * HQ.C;
  else
    HQ.M = -M;
  end if;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, points = {{0, -40}, {50, 40}, {-50, 40}})}));
end Discharge;
