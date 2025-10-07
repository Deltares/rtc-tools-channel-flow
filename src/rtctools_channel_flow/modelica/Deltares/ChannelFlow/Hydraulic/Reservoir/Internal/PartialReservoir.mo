within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;

partial model PartialReservoir
// partial reservoir model for the Hydraulic library: volumes, flows and water levels. 
// The water balance equation is specified here, and the discharge is split into turbine flow and spill flow. 
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  Modelica.Units.SI.Position H;
    // Parameters
  parameter Real Q_nominal=1.0;
equation
  // Water level
  H = HQUp.H;
  // Mass balance
  der(V) / Q_nominal = (HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q)) / Q_nominal;
  // Split outflow between turbine and spill flow
  HQDown.Q + Q_turbine + Q_spill = 0.0;
end PartialReservoir;
