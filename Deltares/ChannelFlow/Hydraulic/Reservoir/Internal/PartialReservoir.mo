within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialReservoir
  extends Deltares.ChannelFlow.Internal.HQCMTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  Modelica.SIunits.Position H;
equation
  // Water level
  H = HQCMUp.H;
  // Mass balance
  der(V) = HQCMUp.Q + HQCMDown.Q + sum(QForcing) + sum(QLateral.Q);
  // Split outflow between turbine and spill flow
  HQCMDown.Q + Q_turbine + Q_spill = 0.0;
end PartialReservoir;
