within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

partial block PartialKNNonlinearVForm
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;

  Modelica.SIunits.Volume V;
  parameter Internal.KNNonlinearityParameterNumerator k_internal_num "Nonlinearity parameter numerator";
  parameter Internal.KNNonlinearityParameterNumerator k_internal_den "Nonlinearity parameter denominator";
  parameter Internal.KNAlpha alpha_internal "Routing parameter";
  parameter SI.Position L;

  parameter Modelica.SIunits.VolumeFlowRate Q0 = 1e-6;

equation
  // We express the storage in terms of the corresponding flows.
  // Note that: V = L * alpha * Q_out ^ k and Q_in - Q_out = der(V).
  // To ensure that V(Q_out) is still differentiable when Q = 0,
  // we approximatetry V with V = L * alpha ((Q_out - Q0) ^ k - Q0 ^ k)) for some small Q0.

  der(V) = QIn.Q - QOut.Q;
  V = L * alpha_internal * (
    (QOut.Q + Q0) ^ (k_internal_num / k_internal_den)
    - Q0 ^ (k_internal_num / k_internal_den)
  );

initial equation
  // Steady state inizialization

  QIn.Q - QOut.Q = 0.0;

end PartialKNNonlinearVForm;
