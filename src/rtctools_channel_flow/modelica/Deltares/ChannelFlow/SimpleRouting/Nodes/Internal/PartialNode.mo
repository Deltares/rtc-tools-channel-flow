within Deltares.ChannelFlow.SimpleRouting.Nodes.Internal;

partial block PartialNode "Partial block with multiple inflows and multiple outflows, where allocation is based on explicitly specified outflows."
  import SI = Modelica.SIunits;
  replaceable parameter Integer nout(min = 0) = 0 "Number of outflows";
  parameter Integer nin(min = 1) = 1 "Number of inflows.";
  Deltares.ChannelFlow.Interfaces.QInPort QIn[nin] annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Deltares.ChannelFlow.Interfaces.QOutPort QOut[nout] annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
protected
  SI.VolumeFlowRate QInSum;
  SI.VolumeFlowRate QOutSum;
equation

  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Text(visible = true, origin = {-80, 40}, extent = {{-20, -20}, {20, 20}}, textString = "%nin"), Text(visible = true, origin = {80, 40}, extent = {{-20, -20}, {20, 20}}, textString = "%nout")}));
end PartialNode;
