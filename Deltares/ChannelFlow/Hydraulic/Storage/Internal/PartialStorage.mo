within Deltares.ChannelFlow.Hydraulic.Storage.Internal;

partial model PartialStorage
  extends Deltares.ChannelFlow.Internal.HQCMOnePort;
  extends Deltares.ChannelFlow.Internal.QForcing;
  extends Deltares.ChannelFlow.Internal.Volume;
  // Homotopy parameter
  parameter Real theta;
  // Nominal values used in linearization
  parameter Modelica.SIunits.Volume V_nominal;
  parameter Modelica.SIunits.Density C_nominal[HQCM.medium.n_substances] = fill(1, HQCM.medium.n_substances);
equation
  der(V) = HQCM.Q + sum(QForcing);
  // der(V*HQ.C) = HQ.M, this equation is linearized.
  HQCM.M = theta * der(V * HQCM.C) + (1 - theta) * (C_nominal * der(V) + V_nominal * der(HQCM.C));
end PartialStorage;
