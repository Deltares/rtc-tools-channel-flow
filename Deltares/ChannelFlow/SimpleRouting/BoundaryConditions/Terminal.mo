within Deltares.ChannelFlow.SimpleRouting.BoundaryConditions;

block Terminal
  output Modelica.SIunits.VolumeFlowRate Q;
  Interfaces.QInPort QIn annotation(Placement(visible = true, transformation(origin = {-80, -0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
equation
  Q = QIn.Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -30}, {50, 30}})}));
end Terminal;