within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir_turbine_out
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q);
  // Outflow is only from the turbine
  QOut.Q = Q_turbine;
end Reservoir_turbine_out;
