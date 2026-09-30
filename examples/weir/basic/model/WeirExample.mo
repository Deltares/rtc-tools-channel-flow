model WeirExample
  input Modelica.SIunits.VolumeFlowRate upstream_q_ext(fixed = true);
  input Modelica.SIunits.VolumeFlowRate downstream_q_ext(fixed = true);
  input Modelica.SIunits.VolumeFlowRate WeirFlow1(fixed = false, nominal=1, min=0, max=2.5);
  output Modelica.SIunits.Position branch_1_water_level;
  output Modelica.SIunits.Position branch_2_water_level;
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge Upstream annotation(Placement(visible = true, transformation(origin = {-90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge Downstream annotation(Placement(visible = true, transformation(origin = {90, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Deltares.ChannelFlow.Hydraulic.Reservoir.Linear Branch1(A = 50, H_b = 0, H(nominal=1, min=0, max=100))  annotation(Placement(visible = true, transformation(origin = {-62, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Deltares.ChannelFlow.Hydraulic.Reservoir.Linear Branch2(A = 100, H_b = 0, H(nominal=1, min=0, max=100))  annotation(Placement(visible = true, transformation(origin = {-8, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Deltares.ChannelFlow.Hydraulic.Structures.Weir weir1(hw_max = 3, hw_min = 1.7, q_max = 1, q_min = 0, width = 10)  annotation(Placement(visible = true, transformation(origin = {-36, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(weir1.HQDown, Branch2.HQUp) annotation(Line(points = {{-28, 10}, {-16, 10}, {-16, 10}, {-16, 10}}, color = {0, 0, 255}));
  connect(Branch1.HQDown, weir1.HQUp) annotation(Line(points = {{-54, 10}, {-44, 10}, {-44, 10}, {-44, 10}}, color = {0, 0, 255}));
  connect(Branch2.HQDown, Downstream.HQ) annotation(Line(points = {{62, 10}, {82, 10}, {82, 10}, {82, 10}}, color = {0, 0, 255}));
  connect(Upstream.HQ, Branch1.HQUp) annotation(Line(points = {{-82, 10}, {-70, 10}, {-70, 10}, {-70, 10}}, color = {0, 0, 255}));
  Branch1.HQDown.H=Branch1.H;
  Branch2.HQDown.H=Branch2.H;
  Branch1.Q_turbine=0;
  Branch2.Q_turbine=0;
  Upstream.Q = upstream_q_ext;
  Downstream.Q = downstream_q_ext;
  WeirFlow1 = weir1.Q;
  branch_1_water_level = Branch1.H;
  branch_2_water_level = Branch2.H;
end WeirExample;
