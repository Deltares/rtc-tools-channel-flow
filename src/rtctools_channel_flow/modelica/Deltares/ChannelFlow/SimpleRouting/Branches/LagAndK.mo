within Deltares.ChannelFlow.SimpleRouting.Branches;

block LagAndK
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  parameter Modelica.Units.SI.Time Lag_parameter = 3600;
  parameter Modelica.Units.SI.Time K_parameter = 1;
  Deltares.ChannelFlow.SimpleRouting.Branches.Delay delay1(duration = Lag_parameter, Q_nominal=Q_nominal) annotation(
    Placement(visible = true, transformation(origin = {-38, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Deltares.ChannelFlow.SimpleRouting.Branches.Muskingum muskingum1(x=0.0, K = K_parameter, Q_nominal=Q_nominal) annotation(
    Placement(visible = true, transformation(origin = {2, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  // Nominal values for scaling
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  QIn.Q / Q_nominal = delay1.QIn.Q / Q_nominal;
  QOut.Q / Q_nominal = muskingum1.QOut.Q / Q_nominal;
  connect(delay1.QOut, muskingum1.QIn) annotation(
    Line(points = {{-30, -6}, {-6, -6}, {-6, -4}, {-6, -4}}));

end LagAndK;
