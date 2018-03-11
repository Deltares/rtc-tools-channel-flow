within Deltares.ChannelFlow.Internal;

partial class QLateral
  parameter Integer n_QLateral(min = 0) = 0;
  Deltares.ChannelFlow.Interfaces.QInPort QLateral[n_QLateral] annotation(Placement(visible = true, transformation(origin = {40, 80}, extent = {{20, -20}, {-20, 20}}, rotation = 90), iconTransformation(origin = {40, 80}, extent = {{20, -20}, {-20, 20}}, rotation = 90)));
  annotation(Icon(graphics = {Text(extent = {{50, 100}, {90, 80}}, textString = "%n_QLateral", fontSize = 80, fontName = "Arial")}, coordinateSystem(initialScale = 0.1)));
end QLateral;
