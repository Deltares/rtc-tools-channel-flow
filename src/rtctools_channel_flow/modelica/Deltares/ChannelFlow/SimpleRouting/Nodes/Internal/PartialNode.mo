within Deltares.ChannelFlow.SimpleRouting.Nodes.Internal;

partial block PartialNode "Partial block with multiple inflows and multiple outflows, where allocation is based on explicitly specified outflows."
  import SI = Modelica.Units.SI;
  replaceable parameter Integer nout(min = 0) = 0 "Number of outflows";
  parameter Integer nin(min = 1) = 1 "Number of inflows.";
  Deltares.ChannelFlow.Interfaces.QInPort QIn[nin](each Q(nominal=Q_nominal)) annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Deltares.ChannelFlow.Interfaces.QOutPort QOut[nout](each Q(nominal=Q_nominal)) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
protected
  SI.VolumeFlowRate QInSum(nominal=Q_nominal);
  SI.VolumeFlowRate QOutSum(nominal=Q_nominal);
equation

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-80, 40}, extent = {{-20, -20}, {20, 20}}, textString = "%nin"), Text(visible = true, origin = {80, 40}, extent = {{-20, -20}, {20, 20}}, textString = "%nout")}));
end PartialNode;
