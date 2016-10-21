within Deltares.ChannelFlow.Internal;

partial class QForcing
  parameter Integer n_QForcing(min = 0) = 0;
  input Modelica.SIunits.VolumeFlowRate QForcing[n_QForcing];
  annotation(Icon(graphics = {Line(origin = {0, 0}, points = {{-30, 90}, {0, 60}, {30, 90}}), Text(origin = {0, 0}, extent = {{-20, 100}, {20, 80}}, textString = "%n_QForcing", fontSize = 80, fontName = "Arial")}));
end QForcing;