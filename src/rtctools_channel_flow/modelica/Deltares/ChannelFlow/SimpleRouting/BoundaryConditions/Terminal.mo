within Deltares.ChannelFlow.SimpleRouting.BoundaryConditions;

block Terminal
  extends Deltares.ChannelFlow.Internal.QSI;
  // Outputs
  output Modelica.SIunits.VolumeFlowRate Q;
equation
  Q = QIn.Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -30}, {50, 30}})}));
end Terminal;
