within Deltares.ChannelFlow.SimpleRouting.Branches;

block Steady
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
equation
  QOut.Q = QIn.Q - QLat.Q + sum(QForcing);
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Line(points = {{-50, 0}, {50, 0}}), Line(origin = {0, 30}, points = {{0, 20}, {0, -30}})}));
end Steady;