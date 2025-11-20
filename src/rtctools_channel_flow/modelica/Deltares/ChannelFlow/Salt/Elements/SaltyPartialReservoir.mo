within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoir
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not user in optimization.
  */
  import SI = Modelica.Units.SI;
  extends QMForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;
  
  SI.Position H;
  SI.Density C(min = 0, nominal = 1);

  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  parameter SI.VolumeFlowRate C_nominal = 1.0;
  parameter SI.Volume V_nominal;

equation
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C[1];
  C = HQDown.C[1];

  // Mass balance, can be used for one medium - salt
  der(V) = HQUp.Q + HQDown.Q + sum(QForcing) + sum(QLateral.Q);
  der(V * C) = HQUp.M[1] + HQDown.M[1] + sum(MForcing[:,1]);
  
end SaltyPartialReservoir;
