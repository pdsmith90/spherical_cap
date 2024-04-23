Pollack 1973 https://doi.org/10.1029/JB078i011p01760

compute gravitational potential of spherical cap for a satellite in 2-body motion, considering only secular j2 perturbation

alt_dep.m plots some simple acceleration vs altitude

spherical_cap_pole.m is a simplified spherical cap located on the pole and thus only contains zonal coefficients.

spherical_cap.m is for arbitrary location cap and field point, includes calculating gravity vector

orbiting_cap.m is wrapper script for spherical_cap.m which can compute potential and gravity anomaly for multiple orbits of a circular satellite, assuming only 2-body and ascending node precession due to J2. 
