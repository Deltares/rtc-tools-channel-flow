within Deltares.ChannelFlow.Salt.Elements;

class QMForcing
  parameter Integer n_QForcing(min = 0) = 0;
  replaceable package medium = Deltares.ChannelFlow.Media.FreshWater;
  input Modelica.SIunits.VolumeFlowRate QForcing[n_QForcing];
  input Modelica.SIunits.VolumeFlowRate MForcing[n_QForcing];
  annotation(Icon(graphics = {Line(origin = {-40, 0}, points = {{-20, 100}, {0, 60}, {20, 100}}), Text(extent = {{-90, 100}, {-50, 80}}, textString = "%n_QForcing")}, coordinateSystem(initialScale = 0.1)));
end QMForcing;
