within Deltares.ChannelFlow.SimpleRouting.Branches;

block Delay
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  parameter SI.Duration duration = 0.0;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  QOut.Q / Q_nominal = delay(QIn.Q, duration) / Q_nominal;
  annotation(Icon(graphics = {Text(extent = {{-25, 25}, {25, -25}}, textString = "τ")}, coordinateSystem(initialScale = 0.1)));
end Delay;