within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialReservoir
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  Modelica.SIunits.Position H;
equation
  // Water level
  H = HQUp.H;
  // Mass balance
  der(V) = HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q);
  // Split outflow between turbine and spill flow
  HQDown.Q + Q_turbine + Q_spill = 0.0;
end PartialReservoir;
