within Deltares.ChannelFlow.Internal;

partial class Reservoir
  import SI = Modelica.SIunits;
  // Inputs
  input SI.VolumeFlowRate Q_turbine;
  input SI.VolumeFlowRate Q_spill;
  // States
  SI.Volume V(min = 0, nominal = 1e6);
equation
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, points = {{40, 50}, {-45, 0}, {40, -50}, {40, 50}, {40, 50}}), Text(origin = {0, -80}, extent = {{-70, 20}, {70, -20}}, textString = "%name", fontName = "MS Shell Dlg 2")}));
end Reservoir;
