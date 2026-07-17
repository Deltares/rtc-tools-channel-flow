within Deltares.ChannelFlow.SimpleRouting.Nodes;

block Node "Block with multiple inflows and multiple outflows and forcing, where allocation is based on explicitly specified outflows."
  import SI = Modelica.Units.SI;
  extends Internal.PartialNode(redeclare parameter Integer nout(min = 1) = 1);
  extends Deltares.ChannelFlow.Internal.QForcing(QForcing(each nominal=Q_nominal));
  input SI.VolumeFlowRate QOut_control[nout - 1](each nominal=Q_nominal);
equation
  QInSum / Q_nominal = sum(QIn.Q) / Q_nominal;
  QOutSum / Q_nominal = (QInSum + sum(QForcing)) / Q_nominal;
  for i in 1:nout - 1 loop
    QOut[i].Q / Q_nominal = QOut_control[i] / Q_nominal;
  end for;
  QOut[nout].Q / Q_nominal = (QOutSum - sum(QOut_control[1:nout - 1])) / Q_nominal;
  annotation(Icon(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10}), graphics = {Polygon(visible = true, fillColor = {255, 170, 0}, fillPattern = FillPattern.Solid, points = {{0, 50}, {-30, 40}, {30, -40}, {0, -50}, {-30, -40}, {30, 40}}), Polygon(visible = true, fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, points = {{-50, 0}, {-40, 30}, {-30, 40}, {30, -40}, {40, -30}, {50, 0}, {40, 30}, {30, 40}, {-30, -40}, {-40, -30}})}), Diagram(coordinateSystem(extent = {{-100, -100}, {100, 100}}, preserveAspectRatio = true, initialScale = 0.1, grid = {10, 10})));
  annotation(
  Documentation(info="
  <html>
  <p>
  A <em>Node</em> connects multiple branches. It can be used to model
  bifurcations and confluences in the water system model.
  </p>
  
  <code>
  Deltares.ChannelFlow.SimpleRouting.Nodes.Node Alder(
    nin = 2,
    nout = 1,
    n_QForcing = 0,
    QIn.Q(each nominal = 0.3),
    QOut.Q(each nominal = 10))
  </code>
  
  <p>
  In this example, the <em>Node</em> with name <em>Alder</em> has two inflow
  connectors and one outflow connector (<code>nout</code>) to represent a confluence. This node is used in an optimization model and nominal values are assigned to the inflow and outflow connectors. It is possible to specify forcing outflow (<code>n_QForcing</code>) to a <em>Node</em>. Typically, forcing is used to model extractions from a channel network. 
  </p>
  

  </html>")
  );
end Node;
