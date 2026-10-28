within Deltares.ChannelFlow.Hydraulic.Reservoir;

model HomotopicVolume
  extends Internal.PartialHomotopicVolume;

  annotation (
    Documentation(info="
  <html>
  <p>
  This model introduces a nonlinear, non-convex relationship between
  reservoir volume and water level using a fourth-order polynomial.
  </p>

  <p>
  The optimization starts from a linear approximation in which the water level
  is computed from the reservoir volume and surface area.
  </p>

  <p>
  During homotopy optimization, the parameter <code>theta</code> gradually
  transitions the model from the linear approximation to the full nonlinear
  polynomial relationship.
  </p>

  <p>
  This model must be used together with the HomotopyMixin.
  </p>
  </html>")
    );

end HomotopicVolume;