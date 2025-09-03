//within Deltares.ChannelFlow.Hydraulic.Structures;
within Deltares.ChannelFlow.Salt.Elements;

model SubstanceControlledStructure "SubstanceControlledStructure"
  extends Deltares.ChannelFlow.Internal.HQTwoPort;
  function smooth_switch = Deltares.ChannelFlow.Internal.Functions.SmoothSwitch;
  function smooth_min = Deltares.ChannelFlow.Internal.Functions.SmoothMin;
  function smooth_abs = Deltares.ChannelFlow.Internal.Functions.SmoothAbs;

  // Homotopy parameter
  parameter Real theta = 1.0;
  // Nominal values used in linearization
  parameter Modelica.SIunits.MassFlowRate Q_nominal = 1;
  parameter Modelica.SIunits.Density C_nominal[HQUp.medium.n_substances] = fill(1e-3, HQUp.medium.n_substances);
  
  SI.Concentration salinity_psu_up(nominal=34.7, start = 34.7);
  SI.Concentration salinity_psu_down(nominal=34.7, start = 34.7);

  parameter SI.Temperature temperature_up;// = 16.0;
  parameter SI.Temperature temperature_down; // = 16.0;

  SI.Density rho_up(nominal=1000, start = 1000.0);
  SI.Density rho_down(nominal=1000, start = 1000.0);

  parameter SI.Height width = 2000;
  
  Real a_up;
  Real b_up;
  Real c_up;
  Real a_down;
  Real b_down;
  Real c_down;
  
  SI.Density rho_ref_up;
  SI.Density rho_ref_down;
  
  Real flux_q1_s1;
  Real epsilon_abs = 0.000001;

equation
  
  salinity_psu_up = HQUp.C[1] / rho_up * 1000.0;
  salinity_psu_down = HQDown.C[1] / rho_down * 1000.0;

  a_up = 8.24493E-1 - 4.0899E-3 * temperature_up + 7.6438E-5 * temperature_up^2.0;// - 8.2467E-7 * temperature_up^3.0 + 5.3875E-9 * temperature_up^4.0;
  b_up = -5.72466E-3 + 1.0227E-4 * temperature_up - 1.6546E-6 * temperature_up^2.0;
  c_up = 4.8314E-4;
  rho_ref_up = (999.842594 + 6.793952E-2 * temperature_up - 9.095290E-3 * temperature_up^2.0 +1.001685E-4 * temperature_up^3.0 - 1.120083E-6 *temperature_up^4.0 +6.536332E-9 * temperature_up^5.0);
  rho_up = rho_ref_up + a_up * salinity_psu_up + b_up * salinity_psu_up^1.5 + c_up * salinity_psu_up^2.0;
  
  a_down = 8.24493E-1 - 4.0899E-3 * temperature_down + 7.6438E-5 * temperature_down^2.0;// - 8.2467E-7 * temperature_down^3.0 + 5.3875E-9 * temperature_down^4.0;
  b_down = -5.72466E-3 + 1.0227E-4 * temperature_down - 1.6546E-6 * temperature_down^2.0;
  c_down = 4.8314E-4;
  rho_ref_down = (999.842594 + 6.793952E-2 * temperature_down - 9.095290E-3 * temperature_down^2.0 +1.001685E-4 * temperature_down^3.0 - 1.120083E-6 *temperature_down^4.0 +6.536332E-9 * temperature_down^5.0); 
  rho_down = rho_ref_down + a_down * salinity_psu_down + b_down * salinity_psu_down^1.5 + c_down * salinity_psu_down^2.0;

  flux_q1_s1 =  (2*9.81)^0.5 * width / 2 * min(HQUp.H, HQDown.H)^1.5*(smooth_abs(rho_up-rho_down, epsilon_abs)/(rho_up+rho_down))^0.5;

  /*
  if HQUp.Q  > flux_q1_s1 then
      HQUp.M = HQUp.Q * HQUp.C[1];
  else
      HQUp.M =0.5 * HQUp.Q * (HQUp.C[1]+ HQDown.C[1]) + (HQUp.C[1]-HQDown.C[1])* 0.5 * (2*9.81)^0.5 * width / 2 * min(HQUp.H, HQDown.H)^1.5*(smooth_abs(rho_up-rho_down, epsilon_abs)/(rho_up+rho_down))^0.5;
  end if;
  */
  
  
  if HQUp.Q < -flux_q1_s1 then
     HQUp.M = HQUp.Q * HQDown.C[1];
  elseif HQUp.Q  > flux_q1_s1 then
      HQUp.M = HQUp.Q * HQUp.C[1];
  else
      HQUp.M =0.5 * HQUp.Q * (HQUp.C[1]+ HQDown.C[1]) + (HQUp.C[1]-HQDown.C[1])* 0.5 * (2*9.81)^0.5 * width / 2 * min(HQUp.H, HQDown.H)^1.5*(smooth_abs(rho_up-rho_down, epsilon_abs)/(rho_up+rho_down))^0.5;
  end if;

 
  annotation(Icon(coordinateSystem( initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(origin = {0, -16.67}, fillColor = {255, 128, 0}, fillPattern = FillPattern.Solid, points = {{0, 66.667}, {-50, -33.333}, {50, -33.333}, {0, 66.667}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
end SubstanceControlledStructure;