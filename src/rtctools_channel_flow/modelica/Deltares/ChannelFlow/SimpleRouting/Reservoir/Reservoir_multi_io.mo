within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir_multi_io "Reservoir with multiple inflows/outflows."
  import SI = Modelica.Units.SI;
  // Parameters
  parameter Integer nin(min = 0) = 1 "Number of inflows.";
  parameter Integer nout(min = 0) = 1 "Number of outflows";
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  // Inflow and outflow
  Deltares.ChannelFlow.Interfaces.QInPort QIn[nin](each Q(nominal=Q_nominal)) annotation(Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Deltares.ChannelFlow.Interfaces.QOutPort QOut[nout](each Q(nominal=Q_nominal)) annotation(Placement(visible = true, transformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  input SI.VolumeFlowRate QOut_control[nout](each nominal=Q_nominal);
  // States
  SI.Volume V(min = 0, nominal = 1e6);
equation
  // QOut from control input.
  for i in 1:nout loop
    QOut[i].Q / Q_nominal = QOut_control[i] / Q_nominal;
  end for;
  // Mass balance
  der(V) / Q_nominal = (sum(QIn.Q) - sum(QOut.Q)) / Q_nominal;  
  // Annotation
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(fillColor = {0, 255, 255}, fillPattern = FillPattern.Solid, points = {{40, 50}, {-45, 0}, {40, -50}, {40, 50}, {40, 50}}), Text(origin = {0, -80}, extent = {{-70, 20}, {70, -20}}, textString = "%name", fontName = "MS Shell Dlg 2")}));
end Reservoir_multi_io;
