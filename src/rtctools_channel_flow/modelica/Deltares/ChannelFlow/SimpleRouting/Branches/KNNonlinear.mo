within Deltares.ChannelFlow.SimpleRouting.Branches;

block KNNonlinear "K-N non-inear routing"
  import SI = Modelica.Units.SI;
  extends Internal.PartialKNNonlinear(k_internal_num=k_num, k_internal_den=k_den, alpha_internal=alpha, L_internal=L);
  parameter Internal.KNNonlinearityParameterNumerator k_num "Nonlinearity parameter numerator";
  parameter Internal.KNNonlinearityParameterDenominator k_den "Nonlinearity parameter denominator";
  parameter Internal.KNAlpha alpha "Routing parameter";
  parameter SI.Position L;
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Rectangle( extent={{-10,10},{10,-10}}, lineColor={0,0,255}, fillColor={85,170,255}, fillPattern=FillPattern.Solid), Line(points = {{-50, 0}, {50, 0}})}));
end KNNonlinear;
