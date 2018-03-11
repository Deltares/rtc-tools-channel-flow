within Deltares.ChannelFlow.SimpleRouting.Reservoir;

block Reservoir
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
equation
  // Mass balance
  der(V) = QIn.Q - QOut.Q + sum(QForcing) + sum(QLateral.Q);
  // Split outflow between turbine and spill flow
  QOut.Q = Q_turbine + Q_spill;
end Reservoir;
