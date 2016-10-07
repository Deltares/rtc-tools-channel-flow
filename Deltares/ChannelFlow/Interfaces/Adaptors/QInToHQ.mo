within Deltares.ChannelFlow.Interfaces.Adaptors;

model QInToHQ "Model with QOutPort to model with HQPort"
  QInPort QIn annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  HQPort HQ annotation(Placement(visible = true, transformation(origin = {80, -0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  output Modelica.SIunits.Position H "Level above datum";
equation
  HQ.Q + QIn.Q = 0;
  H = HQ.H;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, points = {{-50, 0}, {50, 0}})}));
end QInToHQ;