within Deltares.ChannelFlow.Hydraulic.Storage.Internal;

partial model PartialStorage
  extends Deltares.ChannelFlow.Internal.HQOnePort(HQ.C(nominal = C_nominal));
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.Volume;
  // Homotopy parameter
  parameter Real theta = 1.0;
  // Nominal values used in linearization
  parameter Modelica.SIunits.Volume V_nominal;
  parameter Modelica.SIunits.Density C_nominal[HQ.medium.n_substances] = fill(1e-3, HQ.medium.n_substances);
equation
  der(V) = HQ.Q + sum(QForcing);
  // der(V*HQ.C) = HQ.M, this equation is linearized.
  HQ.M = theta * der(V * HQ.C) + (1 - theta) * (C_nominal * der(V) + V_nominal * der(HQ.C));
end PartialStorage;
