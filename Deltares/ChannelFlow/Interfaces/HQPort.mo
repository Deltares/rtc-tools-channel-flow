within Deltares.ChannelFlow.Interfaces;

connector HQPort "Connector with potential water level (H) and flow discharge (Q)"
  Modelica.SIunits.Position H "Level above datum";
  flow Modelica.SIunits.VolumeFlowRate Q "Volume flow (positive inwards)";
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}), Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}})}));
end HQPort;