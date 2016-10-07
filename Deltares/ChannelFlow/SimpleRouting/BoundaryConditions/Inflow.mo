within Deltares.ChannelFlow.SimpleRouting.BoundaryConditions;

block Inflow
  input Modelica.SIunits.VolumeFlowRate Q;
  Interfaces.QOutPort QOut annotation(Placement(visible = true, transformation(origin = {80, -0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
equation
  QOut.Q = Q;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-15, 15}, {-50, 0}, {-15, -15}, {0, -50}, {15, -15}, {50, 0}, {15, 15}})}));
end Inflow;