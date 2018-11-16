within Deltares;

package Constants
  "Library of mathematical constants and constants of nature"

  // Mathematical constants
  final constant Real pi=3.14159265358979;
  final constant Real D2R=0.01745329251994329577 "Degree to Radian";
  final constant Real R2D=57.2957795130823208768 "Radian to Degree";

  // Constants of nature
  // (name, value, description from http://physics.nist.gov/cuu/Constants/index.html, data from 2014)
  final constant Modelica.SIunits.Acceleration g_n=9.80665;

  // Numerical tuning constants
  constant Real eps=1e-12 "Small number used to guard against singularities";
end Constants;