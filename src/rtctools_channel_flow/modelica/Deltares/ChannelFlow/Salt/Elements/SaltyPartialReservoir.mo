//within Deltares.ChannelFlow.Hydraulic.Reservoir.Internal;
within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoir
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends QMForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  Modelica.SIunits.Position H;
  SI.Concentration C(nominal = 0.001);

  parameter Real theta = 1.0;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  parameter SI.VolumeFlowRate C_nominal = 1.0;
  parameter Modelica.SIunits.Volume V_nominal;

equation
  // Water level
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C;
  C = HQDown.C;

  // Mass balance
  der(V) = HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q);
  der(V * C) = HQUp.M + HQDown.M + sum(MForcing);
  //HQ.M / (Q_nominal * C_nominal) = (theta * der(V * HQ.C) + (1 - theta) * Q_nominal * der(HQ.C))  / (Q_nominal * C_nominal);

  // Split outflow between turbine and spill flow
  // HQDown.Q + Q_turbine + Q_spill = 0.0;
end SaltyPartialReservoir;
