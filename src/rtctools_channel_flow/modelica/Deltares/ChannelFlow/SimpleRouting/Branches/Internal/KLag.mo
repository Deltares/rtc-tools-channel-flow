within Deltares.ChannelFlow.SimpleRouting.Branches.Internal;

block KLag "K-lag routing"
  extends PartialKLag(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L=L);
  parameter KLagNonlinearityParameterNumerator k_num = 1 "Nonlinearity parameter numerator";
  parameter KLagNonlinearityParameterDenominator k_den = 1 "Nonlinearity parameter denominator";
  parameter KLagAlpha alpha = 36 "Routing parameter";
  parameter SI.Position L = 100;
end KLag;
