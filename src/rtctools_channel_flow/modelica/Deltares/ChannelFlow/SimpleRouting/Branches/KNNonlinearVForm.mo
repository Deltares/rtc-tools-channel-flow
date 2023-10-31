within Deltares.ChannelFlow.SimpleRouting.Branches;

// K-N non-linear routing using a volume formulation.
// The KNNonlinearVForm formulation seems more stable than that of KNNonlinear.
block KNNonlinearVForm
  import SI = Modelica.SIunits;
  extends Internal.PartialKNNonlinearVForm(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L=L);
  parameter Internal.KNNonlinearityParameterNumerator k_num "Nonlinearity parameter numerator";
  parameter Internal.KNNonlinearityParameterDenominator k_den "Nonlinearity parameter denominator";
  parameter Internal.KNAlpha alpha "Routing parameter";
  parameter SI.Position L;
end KNNonlinearVForm;
