within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QLateral(QLateral.Q(each nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.Reservoir(Q_turbine(nominal=Q_nominal), Q_spill(nominal=Q_nominal));
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // Mass balance
  der(V) / Q_nominal = (QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;
  // Split outflow between turbine and spill flow
  QOut.Q / Q_nominal = (Q_turbine + Q_spill) / Q_nominal;
end Reservoir;
