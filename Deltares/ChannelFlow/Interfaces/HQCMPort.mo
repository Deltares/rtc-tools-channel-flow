within Deltares.ChannelFlow.Interfaces;

connector HQCMPort
  extends HQPort;
  replaceable package medium = Deltares.ChannelFlow.Media.FreshWater;  
  flow Modelica.SIunits.MassFlowRate M[medium.n_substances](each nominal = 10);
  Modelica.SIunits.Density C[medium.n_substances](each min = 0, each nominal = 1);
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {255, 255, 255}, extent = {{-100, -100}, {100, 100}}), Ellipse(visible = true, lineColor = {0, 0, 255}, fillColor = {0, 0, 255}, fillPattern = FillPattern.Solid, extent = {{-50, -50}, {50, 50}})}));
end HQCMPort;