model Example
  import SI = Modelica.Units.SI;
  model MyPumpingStation
    extends Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.PumpingStation(
      n_pumps=1
    );

    Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.Pump pump1(
      power_coefficients = {{{    3934.844148583077,      11874.52561417293,       5631.284774553001},
                      {    5912.865469672678,      -7108.321077407867,            0},
                      {    2243.19274208466,            0,            0}}},

      working_area = {{{         0,            0,            0},
                 {         -1,            0,            0},
                 {         0,            0,            0}},
                {{       0.63713,            -1,            0},
                 {         0,            0,            0},
                 {         0,            0,            0}},
                {{         0,            -1,            0},
                 {         0,            0,            0},
                 {         0,            0,            0}},
                {{      -2.83029,          1.150001,            0},
                 {         -1,            0,            0},
                 {         0,            0,            0}},
                {{         -2,            0,            0},
                 {         -1,            0,            0},
                 {         0,            0,            0}}},


      working_area_direction = {1, 1, -1, -1, -1},
      speed_coefficients = {{      -0.01208409591200699,         59.21659724237954,         -0.0004409413982930602},
                      {      -0.01674254053715328,         -0.007497393806440563,            0},
                      {      -0.005545864740397997,            0,            0}},
      head_option = -1
    ) annotation(Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.Resistance resistance1(C=1.0) annotation(Placement(visible = true, transformation(origin = {-30, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(HQUp, resistance1.HQUp) annotation(Line(points = {{-80, 0}, {-36, 0}, {-36, 0}, {-38, 0}}, color = {0, 0, 255}));
    connect(resistance1.HQDown, pump1.HQUp) annotation(Line(points = {{-22, 0}, {18, 0}, {18, 0}, {18, 0}}, color = {0, 0, 255}));
    connect(pump1.HQDown, HQDown) annotation(Line(points = {{34, 0}, {80, 0}}, color = {0, 0, 255}));
  end MyPumpingStation;

 // Elements in model flow chart
  Deltares.ChannelFlow.Hydraulic.Storage.Linear storage(
    A = 149000,
    H_b = -1.0,
    HQ.H(min = -0.7, max = 0.2),
    V(nominal = 1E5)
  ) annotation(Placement(visible = true, transformation(origin = {-28, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Level sea annotation(Placement(visible = true, transformation(origin = {66, 10}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Deltares.ChannelFlow.Hydraulic.BoundaryConditions.Discharge inflow annotation(Placement(visible = true, transformation(origin = {-72, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  MyPumpingStation pumpingstation1 annotation(Placement(visible = true, transformation(origin = {20, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  // Input variables
  input SI.VolumeFlowRate Q_in(fixed = true);
  input SI.Position H_ext(fixed=true);

  // Energy price is typically of units EUR/kWh (when optimizing for energy
  // usage), but one can also choose for e.g. ton CO2/kWh to get the lowest
  // CO2 output.
  input Real energy_price(fixed=true);

  // NOTE: Because we cannot flag each pump's .Q as "input", we need an extra
  // variable to do this. Format is expected to be the fully specified name,
  // with all dots replaced with underscores.
  input Real pumpingstation1_pump1_Q;
  // TODO: Move bounds to the mixin.
  input Real pumpingstation1_resistance1_dH(min=0.0, max=10.0);

  // Output variables
  output SI.Position storage_level;
  output SI.Position sea_level;
equation
  connect(pumpingstation1.HQUp, storage.HQ) annotation(Line(points = {{12, 10}, {-20, 10}}, color = {0, 0, 255}));
  connect(pumpingstation1.HQDown, sea.HQ) annotation(Line(points = {{28, 10}, {58, 10}}, color = {0, 0, 255}));
  connect(inflow.HQ, storage.HQ) annotation(Line(points = {{-64, 10}, {-20, 10}}, color = {0, 0, 255}));
  // Mapping of variables
  inflow.Q = Q_in;
  sea.H = H_ext;
  pumpingstation1.pump1.Q = pumpingstation1_pump1_Q;
  pumpingstation1.resistance1.dH = pumpingstation1_resistance1_dH;
  storage_level = storage.HQ.H;
  sea_level = H_ext;
end Example;
