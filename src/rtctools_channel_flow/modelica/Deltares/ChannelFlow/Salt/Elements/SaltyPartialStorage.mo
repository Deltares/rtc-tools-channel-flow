within Deltares.ChannelFlow.Salt.Elements;

partial model SaltyPartialStorage
  /*
  This block is designed to be used together with the "salt_simulation_mixin" to calculate dispersive and advective transport
   between salty reservoir elements, do not user in optimization.
  */
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.HQOnePort(HQ.Q(each nominal = Q_nominal), HQ.M(each nominal =Q_nominal - C_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal = Q_nominal));
  extends Deltares.ChannelFlow.Internal.Volume;

  parameter SI.Volume V_nominal;
  parameter SI.Density C_nominal = 1e-3;
  parameter SI.VolumeFlowRate Q_nominal = 1.0;

equation
  
  der(V) / Q_nominal = (HQ.Q + sum(QForcing)) / Q_nominal;
  HQ.M / (Q_nominal * C_nominal) =  der(V * HQ.C)   / (Q_nominal * C_nominal);
  
end SaltyPartialStorage;
