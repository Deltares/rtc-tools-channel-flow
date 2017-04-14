within Deltares.ChannelFlow.Internal;

partial model HQCMTwoPort
  replaceable package medium = Deltares.media.FreshWater;
  Deltares.ChannelFlow.Interfaces.HQCMPort HQCMUp(redeclare package medium = medium) annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Deltares.ChannelFlow.Interfaces.HQCMPort HQCMDown(redeclare package medium = medium) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)), Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {2, 2})));
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Line(visible = true, origin = {-80, -0}, points = {{-20, 30}, {20, 0}, {-20, -30}}, color = {0, 0, 255}), Line(visible = true, origin = {80, 0}, points = {{-20, 30}, {20, 0}, {-20, -30}}, color = {0, 0, 255})}));
end HQCMTwoPort;
