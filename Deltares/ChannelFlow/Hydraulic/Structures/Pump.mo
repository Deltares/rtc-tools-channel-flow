within Deltares.ChannelFlow.Hydraulic.Structures;

model Pump "Pump"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  input Modelica.SIunits.VolumeFlowRate Q;
equation
  HQUp.Q + HQDown.Q = 0;
  HQUp.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, origin = {0, -16.667}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, lineThickness = 2, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end Pump;
