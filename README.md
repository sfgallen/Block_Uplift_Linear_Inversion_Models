# Linear Inversion of River Profiles Assuming Block Uplift

This repository contains two codes that preform different versions of the linear river profile inversion method described in Goren et al. (2014) for spatially uniform, temporally varying (block) uplift assuming uniform erodibility. The method relies on the detachment-limited stream power equation (e.g. Whipple and Tucker, 1999) that describes incision into bedrock, E, as power functions of local channel slope, S, and upstream contributing drainage area, A as:

E = K*A^m*S^n

where m and n are positive constants.

Both sets of codes assume that n = 1 and the m/n ratio is 0.5 unless otherwise defined by the user.

The codes in the “nondimensional_block_uplift” folder follow Goren et al.’s approach to using the transport distance variable chi and river network elevations to derive a “nondimensional uplift” term (note if a reference drainage area of 1 m^2 is used, this is the average normalized steepness index (ksn) per unit chi).

The codes in the “dimensional_block_uplift” folder require an input erodibility, K, parameter to convert chi to the river response time, tau and derive uplift rates in natural units.

Both master functions have help menus, and users are directed to Goren et al. (2014) and Gallen (2018) for details on the methodology and example applications.

Contact:
sean.gallen[at]colostate[dot]edu

# References:

Gallen, S.F., 2018. Lithologic controls on landscape dynamics and aquatic species evolution in post-orogenic mountains. Earth and Planetary Science Letters, 493, pp.150-160.

Goren, L., Fox, M. and Willett, S.D., 2014. Tectonics from fluvial topography using formal linear inversion: Theory and applications to the Inyo Mountains, California. Journal of Geophysical Research: Earth Surface, 119(8), pp.1651-1681.

Whipple, K.X. and Tucker, G.E., 1999. Dynamics of the stream‐power river incision model: Implications for height limits of mountain ranges, landscape response timescales, and research needs. Journal of Geophysical Research: Solid Earth, 104(B8), pp.17661-17674.
