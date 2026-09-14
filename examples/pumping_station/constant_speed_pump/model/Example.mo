model Example
  import SI = Modelica.Units.SI;
  model MyPumpingStation
    extends Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.PumpingStation(
      n_pumps=1
    );

    Deltares.ChannelFlow.Hydraulic.Structures.PumpingStation.Pump pump1(
      power_coefficients = {{{-39095.7484857 ,  66783.26438029},
                             {  7281.58399033,       0.0}},
                            {{-27859.72396609,  52991.89678362},
                             {  6547.78028277,       0.0}},
                            {{-23247.11759041,  47371.07046571},
                             {  6106.48287875,       0.0}},
                            {{-33451.06236192,  59110.131313  },
                             {  7034.0025483,       0.0}},
                            {{-39861.93265223,  67925.88638803},
                             {  7281.81109575,       0.0}},
                            {{-18905.97159263,  40725.7120928 },
                             {  5785.42670305,       0.0}},
                            {{-69080.26073112,  97844.92974535},
                             {  8776.27578528,       0.0}},
                            {{-27357.47031974,  52938.41395038},
                             {  6406.26137405,       0.0}},
                            {{-49670.43493511,  78653.65042992},
                             {  7840.68764271,       0.0}},
                            {{ -8104.87351317,  21084.47678188},
                             {  4586.35220572,       0.0}},
                            {{-35797.98792283,  54481.08143225},
                             {  7595.20055527,       0.0}},
                            {{-60204.93398172,  89038.39027083},
                             {  8408.00958411,       0.0}},
                            {{-22131.50613895,  39713.73557364},
                             {  6465.17137718,       0.0}},
                            {{-14217.62875014,  34100.99360442},
                             {  5201.88819881,       0.0}},
                            {{-65348.50815507,  94475.38978919},
                             {  8575.77174818,       0.0}},
                            {{-39565.60866869,  63810.57380161},
                             {  7597.68369003,       0.0}},
                            {{-15114.73306987,  33754.82613809},
                             {  5508.58396756,       0.0}},
                            {{-51205.68182194,  79946.23472306},
                             {  7963.3995412,       0.0}},
                            {{ -6701.19277067,  20672.59101587},
                             {  4071.93143896,       0.0}},
                            {{-30637.27553712,  53847.70016923},
                             {  6976.92174728,       0.0}},
                            {{-23174.20780931,  45291.10739497},
                             {  6315.4012685,       0.0}},
                            {{-16838.11645224,  33397.62607586},
                             {  5890.34692945,       0.0}},
                            {{-48684.40184549,  73130.85845968},
                             {  8148.05117355,       0.0}},
                            {{ -5585.5596654 ,  17917.6377768 },
                             {  3903.50760056,       0.0}},
                            {{-11530.80023073,  25643.30751799},
                             {  5209.61580979,       0.0}},
                            {{-54050.93126723,  81474.03659672},
                             {  8244.19508265,       0.0}},
                            {{-17059.63308541,  38705.96069607},
                             {  5468.15043893,       0.0}},
                            {{-42990.53693314,  61384.80799185},
                             {  8077.30306983,       0.0}},
                            {{-28740.24310978,  47105.49859182},
                             {  7058.4428374,       0.0}},
                            {{ -9433.46984887,  26148.5409928 },
                             {  4488.53180172,       0.0}},
                            {{-54478.410763  ,  83635.76502236},
                             {  8076.34522102,       0.0}},
                            {{-34733.3023365 ,  61966.5148778 },
                             {  6963.68444901,       0.0}},
                            {{ -4371.73397616,  13785.26261371},
                             {  3759.68617973,       0.0}},
                            {{-44469.5263585 ,  71683.6367841 },
                             {  7713.90758475,       0.0}},
                            {{-10835.36867811,  27900.70743971},
                             {  4844.53856386,       0.0}}},

      working_area={{{  16.6564650552,    -13.8900458136},
                     {  -1.0,               0.0}},
                    {{  16.6564650552,    -13.8900458136},
                     {  -1.0,               0.0}},
                    {{  -9.77101908988,    40.0404314726},
                     {  -1.0,               0.0}},
                    {{  -0.75242261,        2.79784166},
                     {  -1.0,               0.0}},
                    {{  -0.50567904,        2.3278547},
                     {  -1.0,               0.0}},
                    {{  -1.31944138,        3.6515112},
                     {  -1.0,               0.0}},
                    {{  -1.93409137,        4.35798345},
                     {  -1.0,               0.0}},
                    {{  -1.04947144,        3.27458583},
                     {  -1.0,               0.0}}},

      speed_coefficients = {{ 585.0 }},

      working_area_direction = {-1, 1, 1, -1, -1, -1, -1, -1},

      minimum_on=3.0*3600
    ) annotation(Placement(visible = true, transformation(origin = {26, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(HQUp, pump1.HQUp) annotation(Line(points = {{-80, 0}, {18, 0}}, color = {0, 0, 255}));
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
  storage_level = storage.HQ.H;
  sea_level = H_ext;
end Example;
