# Linear Inversion of River Profiles Assuming Block Uplift

## Author: Sean F. Gallen, 2020. [![DOI](https://zenodo.org/badge/314364091.svg)](https://zenodo.org/badge/latestdoi/314364091)

This repository contains two codes that preform different versions of the linear river profile inversion method described in Goren et al. (2014) for spatially uniform, temporally varying (block) uplift assuming uniform erodibility. The method relies on the detachment-limited stream power equation (e.g. Whipple and Tucker, 1999) that describes incision into bedrock, E, as power functions of local channel slope, S, and upstream contributing drainage area, A as:

E = K*A^m*S^n

where m and n are positive constants.

The codes rely on Matlab-based TopoToolbox (Schwanghart and Scherler, 2014; https://topotoolbox.wordpress.com/) to analyze a digital elevation model and define the stream network.

Both sets of codes assume that n = 1 and the m/n ratio is 0.5 unless otherwise defined by the user.

The codes in the “nondimensional_block_uplift” folder follow Goren et al.’s approach to using the transport distance variable chi and river network elevations to derive a “nondimensional uplift” term (note if a reference drainage area of 1 m^2 is used, this is the average normalized steepness index (ksn) per unit chi).

The codes in the “dimensional_block_uplift” folder require an input erodibility, K, parameter to convert chi to the river response time, tau and derive uplift rates in natural units.

Both master functions have help menus, and users are directed to Goren et al. (2014), Gallen (2018), and Pavano and Gallen (in review) for details on the methodology and example applications. Please cite the above mentioned studies if these codes are used for research.

Contact:
sean.gallen[at]colostate[dot]edu

# References:

Gallen, S.F., 2018. Lithologic controls on landscape dynamics and aquatic species evolution in post-orogenic mountains. Earth and Planetary Science Letters, 493, pp.150-160.

Goren, L., Fox, M. and Willett, S.D., 2014. Tectonics from fluvial topography using formal linear inversion: Theory and applications to the Inyo Mountains, California. Journal of Geophysical Research: Earth Surface, 119(8), pp.1651-1681.

Pavano and Gallen, in review. A Geomorphic Examination of the Calabrian Forearc Translation: for consideration in Tectonics.

Schwanghart, W. and Scherler, D., 2014. TopoToolbox 2–MATLAB-based software for topographic analysis and modeling in Earth surface sciences. Earth Surface Dynamics, 2(1), pp.1-7.

Whipple, K.X. and Tucker, G.E., 1999. Dynamics of the stream‐power river incision model: Implications for height limits of mountain ranges, landscape response timescales, and research needs. Journal of Geophysical Research: Solid Earth, 104(B8), pp.17661-17674.
