within Deltares.ChannelFlow.SimpleRouting.Branches;

block KNNonlinear "K-N non-inear routing"
  import SI = Modelica.SIunits;
  extends Internal.PartialKNNonlinear(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L=L);
  parameter Internal.KNNonlinearityParameterNumerator k_num "Nonlinearity parameter numerator";
  parameter Internal.KNNonlinearityParameterDenominator k_den "Nonlinearity parameter denominator";
  parameter Internal.KNAlpha alpha "Routing parameter";
  parameter SI.Position L;
end KNNonlinear;
