within Deltares.ChannelFlow.Interfaces;

connector QInPort "Connector with potential water level (H) and flow discharge (Q)"
  input Modelica.SIunits.VolumeFlowRate Q "Volume flow";
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {0, 255, 0}, fillPattern = FillPattern.Solid, points = {{-100, 100}, {-100, -100}, {100, 0}, {-100, 100}})}));
end QInPort;