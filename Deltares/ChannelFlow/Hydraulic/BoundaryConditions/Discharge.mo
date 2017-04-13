within Deltares.ChannelFlow.Hydraulic.BoundaryConditions;

model Discharge "Defines a discharge"
  extends Deltares.ChannelFlow.Internal.HQCMOnePort;
  input Modelica.SIunits.VolumeFlowRate Q;
  input Modelica.SIunits.MassFlowRate M[HQCM.Mediumport.n_substances];
  parameter Boolean upwind = true; //if true and there is outlfow from the system (inot the discharge boudnary) then the concentration of the connectpoint (system) is used.
equation
  HQCM.Q + Q = 0;
  if upwind then
	if Q > 0 then
		HQCM.M = -M;
	else
		HQCM.M = -Q * HQCM.C;
	end if;
  else
	HQCM.M = - M;
  end if;
  
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 0, 255}, fillPattern = FillPattern.Solid, points = {{0, -40}, {50, 40}, {-50, 40}})}));
end Discharge;
