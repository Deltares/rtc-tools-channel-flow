within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoirBnd
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  SI.Position H;
  SI.Concentration C(nominal = 0.001);

  parameter Real theta = 1.0;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  parameter SI.VolumeFlowRate C_nominal = 1.0;
  parameter SI.Volume V_nominal;

equation
  // Water level
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C[1];
  C = HQDown.C[1];
  //HQUp.Q + HQDown.Q = 0;
  //HQUp.M + HQDown.M = 0;
  // Mass balance
  der(V) = 0; // + sum(QForcing) + sum(QLateral.Q);
  der(V * C) = 0;
  //HQ.M / (Q_nominal * C_nominal) = (theta * der(V * HQ.C) + (1 - theta) * Q_nominal * der(HQ.C))  / (Q_nominal * C_nominal);

  // Split outflow between turbine and spill flow
  // HQDown.Q + Q_turbine + Q_spill = 0.0;
end SaltyPartialReservoirBnd;
