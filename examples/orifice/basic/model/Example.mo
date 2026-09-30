model MyPumpingStation
  extends Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.PumpingStation(
    n_pumps=1
  );

  Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.Pump pump1(
    power_coefficients = {{{3522.8, -27946.3, 54484.8},
                          {1665.43, 5827.81,     0.0},
                          {208.251,     0.0,     0.0}}},

    working_area = {{{ -5.326999,  54.050758,  0.000000},
                     { -1.0,        0.0,       0.0}},
                    {{  0.000426,  -0.001241,  2.564056},
                     { -1.0,        0.0,       0.0}},
                    {{  2.577975,  -5.203480,  0.000000},
                     { -1.0,        0.0,       0.0}},
                    {{ 13.219650,  -3.097600, -7.551339},
                     { -1.0,        0.0,       0.0}}},

    speed_coefficients = {{88.3224,  281.642,  143.011},
                          {58.8183,  -24.9845,   0.0},
                          {-1.45483,   0.0,      0.0}},

    working_area_direction = {1, -1, -1, 1},

    minimum_on=0
  ) annotation(Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(HQUp, pump1.HQUp) annotation(Line(points = {{-80, 0}, {-36, 0}, {-36, 0}, {-38, 0}}, color = {0, 0, 255}));
    connect(pump1.HQDown, HQDown) annotation(Line(points = {{34, 0}, {80, 0}}, color = {0, 0, 255}));
end MyPumpingStation;


model Example
  // Elements in model flow chart
  Deltares.ChannelFlow.Hydraulic.Storage.Linear storage(
    A = 149000,
    H_b = -1.0,
    HQ.H(min=-0.5, max=0.2),
    V(nominal=1E5)
  ) annotation(Placement(visible = true, transformation(origin = {-28, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level sea annotation(Placement(visible = true, transformation(origin = {66, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge inflow annotation(Placement(visible = true, transformation(origin = {-72, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Deltares.ChannelFlow.Hydraulic.Structures.Orifice orifice1(dH_max=2, area=2.4) annotation(Placement(visible = true, transformation(origin = {-28, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));

  MyPumpingStation pumpingstation1 annotation(Placement(visible = true, transformation(origin = {20, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  // Input variables
  input Modelica.SIunits.VolumeFlowRate Q_in(fixed = true);
  input Modelica.SIunits.Position H_ext(fixed=true);
  input Real energy_price(fixed=true);

  // Control variables
  // NOTE: Because we cannot flag the pump's .Q or the orifice's Q as "input",
  // we need extra variables to do this. Format is expected to be the fully
  // specified name, with all dots replaced with underscores.
  input Modelica.SIunits.VolumeFlowRate pumpingstation1_pump1_Q;
  input Modelica.SIunits.VolumeFlowRate orifice1_Q = orifice1.Q;

  // Output variables
  output Modelica.SIunits.Position storage_level;
  output Modelica.SIunits.Position sea_level;
equation
  connect(pumpingstation1.HQUp, storage.HQ) annotation(Line(points = {{12, 10}, {-20, 10}}, color = {0, 0, 255}));
  connect(pumpingstation1.HQDown, sea.HQ) annotation(Line(points = {{28, 10}, {58, 10}}, color = {0, 0, 255}));
  connect(inflow.HQ, storage.HQ) annotation(Line(points = {{-64, 10}, {-20, 10}}, color = {0, 0, 255}));
  connect(orifice1.HQDown, sea.HQ) annotation(Line(points = {{28, 40}, {43, 40}, {43, 10}, {58, 10}}, color = {0, 0, 255}));
  connect(orifice1.HQUp, storage.HQ) annotation(Line(points = {{12, 40}, {-4, 40}, {-4, 10}, {-20, 10}}, color = {0, 0, 255}));
  // Mapping of variables
  inflow.Q = Q_in;
  sea.H = H_ext;
  pumpingstation1.pump1.Q = pumpingstation1_pump1_Q;
  orifice1.Q = orifice1_Q;
  storage_level = storage.HQ.H;
  sea_level = H_ext;
end Example;
