within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

block KLag "K-lag routing"
  extends PartialKLag(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L=L);
  parameter KLagNonlinearityParameterNumerator k_num "Nonlinearity parameter numerator";
  parameter KLagNonlinearityParameterDenominator k_den "Nonlinearity parameter denominator";
  parameter KLagAlpha alpha "Routing parameter";
  parameter SI.Position L;
end KLag;
