within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialReservoirBnd
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not user in optimization.
  */
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.QLateral;
  extends Deltares.ChannelFlow.Internal.Reservoir;

  SI.Position H;
  SI.Density C[medium.n_substances](each min = 0, each nominal = 1);

  parameter SI.VolumeFlowRate Q_nominal = 1.0;
  parameter SI.VolumeFlowRate C_nominal = 1.0;
  parameter SI.Volume V_nominal;

equation
  
  H = HQUp.H;
  H = HQDown.H;

  C = HQUp.C;
  C = HQDown.C;
  
  der(V) = 0;
  der(V * C) = fill(0.0, medium.n_substances);

end SaltyPartialReservoirBnd;
