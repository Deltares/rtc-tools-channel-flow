within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

partial block PartialKLag
  import SI = Modelica.SIunits;
  extends Deltares.ChannelFlow.Internal.QSISO;

  // Note: correct formulation guaranteed only if implicit_step_size is set to the input step size.
  input SI.Duration implicit_step_size(fixed = true);

  parameter Internal.KLagNonlinearityParameterNumerator k_internal_num "Nonlinearity parameter numerator";
  parameter Internal.KLagNonlinearityParameterNumerator k_internal_den "Nonlinearity parameter denominator";
  parameter Internal.KLagAlpha alpha_internal "Routing parameter";
  parameter SI.Position L;

  input Modelica.SIunits.VolumeFlowRate q_out_prev;
  parameter Real min_divisor = Deltares.Constants.eps;

equation
  // We express the storage in terms of the corresponding flows.
  // Note that: V = L * alpha * Q_out ^ k and Q_in - Q_out = der(V).

// Use same trick as Muskingum

  implicit_step_size * (QIn.Q - QOut.Q) / (L * alpha_internal) = (QOut.Q + min_divisor) ^ (k_internal_num / k_internal_den) - (q_out_prev + min_divisor) ^ (k_internal_num / k_internal_den);

  q_out_prev = QOut.Q - implicit_step_size * der(QOut.Q);


initial equation
  // Steady state inizialization

  QIn.Q - QOut.Q = 0.0;

end PartialKLag;
