within Deltares.ChannelFlow.Internal;

partial class QForcing
  parameter Integer n_QForcing(min = 0) = 0;
  input Modelica.SIunits.VolumeFlowRate QForcing[n_QForcing];
  annotation(Icon(graphics = {Line(origin = {-40, 0}, points = {{-20, 100}, {0, 60}, {20, 100}}), Text(extent = {{-90, 100}, {-50, 80}}, textString = "%n_QForcing", fontSize = 80, fontName = "Arial")}, coordinateSystem(initialScale = 0.1)));
end QForcing;
