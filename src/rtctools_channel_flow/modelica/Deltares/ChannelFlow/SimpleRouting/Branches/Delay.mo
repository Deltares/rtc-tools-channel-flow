within Deltares.ChannelFlow.SimpleRouting.Branches;

block Delay
  extends Deltares.ChannelFlow.Internal.QSISO;
  parameter SI.Duration duration = 0.0;
equation
  QOut.Q = delay(QIn.Q, duration);
  annotation(Icon(graphics = {Text(extent = {{-25, 25}, {25, -25}}, textString = "τ")}, coordinateSystem(initialScale = 0.1)));
end Delay;