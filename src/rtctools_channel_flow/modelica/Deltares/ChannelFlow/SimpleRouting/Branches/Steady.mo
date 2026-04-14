within Deltares.ChannelFlow.SimpleRouting.Branches;

block Steady
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QLateral(QLateral.Q(each nominal=Q_nominal));
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  QOut.Q / Q_nominal = (QIn.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Line(points = {{-50, 0}, {50, 0}})}));
end Steady;
