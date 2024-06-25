within Deltares.ChannelFlow.SimpleRouting.Structures;

block RunOfRiverHydropowerComplexFixedHead "Node for a simple complex of a run-of-river hydropower plant and a weir. Head difference for power production is constant."
  import SI = Modelica.Units.SI;
  extends Deltares.ChannelFlow.Internal.QSISO(QIn.Q(nominal=Q_nominal), QOut.Q(nominal=Q_nominal));
  parameter SI.VolumeFlowRate Q_nominal = 1.0;
    // Head difference
  parameter SI.Position dH;
  // Turbine efficiency
  parameter Real nu;
  // Water density
  parameter SI.Density ro;
  // Turbine flow
  output SI.VolumeFlowRate Q_turbine(min=0);
  // Spill flow
  output SI.VolumeFlowRate Q_spill(min=0);
  // Power production
  output SI.Power P;
  equation
    QOut.Q / Q_nominal = (Q_turbine + Q_spill) / Q_nominal;
    QOut.Q / Q_nominal = QIn.Q / Q_nominal;
    P / (nu * ro * Deltares.Constants.g_n * dH * Q_nominal) = Q_turbine / Q_nominal;
end RunOfRiverHydropowerComplexFixedHead;
