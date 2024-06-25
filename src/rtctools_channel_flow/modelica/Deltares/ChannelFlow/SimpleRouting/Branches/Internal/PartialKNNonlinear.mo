within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

partial block PartialKNNonlinear
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));

  // Note: correct formulation guaranteed only if implicit_step_size is set to the input step size.
  input SI.Duration implicit_step_size(fixed = true);

  parameter Internal.KNNonlinearityParameterNumerator k_internal_num "Nonlinearity parameter numerator";
  parameter Internal.KNNonlinearityParameterNumerator k_internal_den "Nonlinearity parameter denominator";
  parameter Internal.KNAlpha alpha_internal "Routing parameter";
  parameter SI.Position L_internal;

  input Modelica.Units.SI.VolumeFlowRate q_out_prev(nominal=Q_nominal);
  parameter Real min_divisor = Deltares.Constants.eps;

  // Nominal values for scaling
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
equation
  // We express the storage in terms of the corresponding flows.
  // Note that: V = L_internal * alpha * Q_out ^ k and Q_in - Q_out = der(V).

// Use same trick as Muskingum

  implicit_step_size * (QIn.Q - QOut.Q) / (L_internal * alpha_internal) = (QOut.Q + min_divisor) ^ (k_internal_num / k_internal_den) - (q_out_prev + min_divisor) ^ (k_internal_num / k_internal_den);

  q_out_prev / Q_nominal = (QOut.Q - implicit_step_size * der(QOut.Q)) / Q_nominal;


initial equation
  // Steady state inizialization

  QIn.Q - QOut.Q = 0.0;

end PartialKNNonlinear;
