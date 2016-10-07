within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model Discharge "Defines a discharge"
  extends Deltares.ChannelFlow.Internal.HQOnePort;
  input Modelica.SIunits.VolumeFlowRate Q;
equation
  HQ.Q + Q = 0;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, points = {{0, -40}, {50, 40}, {-50, 40}})}));
end Discharge;
