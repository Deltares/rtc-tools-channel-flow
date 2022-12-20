within Deltares.ChannelFlow.SimpleRouting.Branches;

block KLag "K-lag routing"
  import SI = Modelica.SIunits;
  extends Internal.PartialKLag(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L=L);
  parameter Internal.KLagNonlinearityParameterNumerator k_num "Nonlinearity parameter numerator";
  parameter Internal.KLagNonlinearityParameterDenominator k_den "Nonlinearity parameter denominator";
  parameter Internal.KLagAlpha alpha "Routing parameter";
  parameter SI.Position L;
end KLag;
