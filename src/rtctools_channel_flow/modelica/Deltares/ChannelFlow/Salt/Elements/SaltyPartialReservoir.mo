within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoir
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends QMForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  // States
  SI.Position H;
  SI.Concentration C(nominal = 0.001);

  parameter Real theta = 1.0;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  parameter SI.Concentration C_nominal = 1.0;
  parameter SI.Volume V_nominal;

equation
  // Water level
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C[1];
  C = HQDown.C[1];

  // Mass balance
  der(V) = HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q);
  der(V * C) = HQUp.M[1] + HQDown.M[1] + sum(MForcing);

end SaltyPartialReservoir;
