within Deltares.ChannelFlow.SimpleRouting.Branches;

block LagAndK
  extends Deltares.ChannelFlow.Internal.QSISO;
  parameter Modelica.SIunits.Time Lag_parameter = 3600;
  parameter Modelica.SIunits.Time K_parameter = 1;
  Deltares.ChannelFlow.SimpleRouting.Branches.Delay delay1(duration = Lag_parameter) annotation(
    Placement(visible = true, transformation(origin = {-38, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Deltares.ChannelFlow.SimpleRouting.Branches.Muskingum muskingum1(x=0.0, K = K_parameter) annotation(
    Placement(visible = true, transformation(origin = {2, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  QIn.Q = delay1.QIn.Q;
  QOut.Q = muskingum1.QOut.Q;
  connect(delay1.QOut, muskingum1.QIn) annotation(
    Line(points = {{-30, -6}, {-6, -6}, {-6, -4}, {-6, -4}}));

end LagAndK;
