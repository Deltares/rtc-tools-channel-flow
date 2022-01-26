within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir_turbine_out
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.QLateral(QLateral.Q(nominal=Q_nominal));
  extends Deltares.ChannelFlow.Internal.Reservoir(Q_turbine(nominal=Q_nominal), Q_spill(nominal=Q_nominal));
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // Mass balance
  der(V) / Q_nominal = (QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;
  // Outflow is only from the turbine
  QOut.Q = Q_turbine;
end Reservoir_turbine_out;
