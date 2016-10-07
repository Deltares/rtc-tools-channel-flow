within Deltares.ChannelFlow.Interfaces.Adaptors;

model HQToQOut "Model with HQPort to model with QinPort"
  HQPort HQ annotation(Placement(visible = true, transformation(origin = {-80, -0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  QOutPort QOut annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  input Modelica.SIunits.Position H "Level above datum";
equation
  QOut.Q = HQ.Q;
  H = HQ.H;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, points = {{-50, 0}, {50, 0}})}));
end HQToQOut;