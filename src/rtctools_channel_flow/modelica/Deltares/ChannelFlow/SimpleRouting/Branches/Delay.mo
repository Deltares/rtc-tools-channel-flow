within Deltares.ChannelFlow.SimpleRouting.Branches;

block Delay
  extends Deltares.ChannelFlow.Internal.QSISO;
  annotation(Icon(graphics = {Text(extent = {{-25, 25}, {25, -25}}, textString = "τ")}, coordinateSystem(initialScale = 0.1)));
end Delay;