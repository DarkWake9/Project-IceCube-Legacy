
IceCube Upgrade data release v01_00 : neutrino Monte Carlo (MC) data

Questions can be directed to data@icecube.wisc.edu


The IceCube Upgrade:
--------------------

The IceCube Upgrade is an upcoming extension to the IceCube neutrino observatory which will densely 
instrument a 2 Mton region in a central region of the deepest ice (within the existing DeepCore 
sub-array). This new sub-array will comprise 7 new strings featuring multi-PMT optical modules (known 
as mDOMs and DEggs) at a vertical spacing of 3m, in addition to a host of new calibration devices. 

The IceCube Upgrade sub-array will offer significantly enhanced neutrino detection efficiency and 
reconstruction resolution relative to DeepCore, as well as lowering the energy threshold to provide 
sensitivity to O(1 GeV) neutrinos.

See [1] for more details of the IceCube Upgrade detector and its sensitivity to neutrino oscillations.


Contents:
---------

This data release contains MC neutrino events from simulations of the IceCube Upgrade. The event 
generator used is GENIE [2]. The events provided are at analysis level, meaning that background removal 
methods have been utilised to heavily suppress unwanted atmospheric muon and coincident detector noise 
trigger events. This data release only includes events from within a 2 Mton IceCube Upgrade cylindrical 
fiducial region enclosing the new sensors.

The data is provided in two alternative formats (both contain the same data, e.g. one or the other 
should be used):
- Individual events
- Binned event distributions, including effective area

Note that data provided in the "individual events" is provided in a similar (but not identical) format 
to the "Three-year high-statistics neutrino oscillation samples" data release using IceCube/DeepCore 
data/simulations (https://icecube.wisc.edu/science/data/highstats_nuosc_3y).

Files:

- README.txt        : This document
- geom_surface.pdf  : Image showing the positions of the IceCube Upgrade strings from a plan/surface
                      view. Nearby IceCube/DeepCore strings are also shown.
- geom_depth.pdf    : Image showing the depths of the sensors and other devices on the IceCube Upgrade 
                      strings. The region shown is the densely instrumented "physics region" in the 
                      deep ice. An IceCube and a DeepCore string are also shown for reference. 
                      Note that the numbers of the x axis are string ID numbers.

- events
-- neutrino_mc.csv   : CSV file containing individual simulated neutrino events. A single file contains 
                       all events.
-- example.ipynb     : A notebook giving examples of how to use this these events

- binned            
-- nu*.csv           : CSV files containing binned simulated neutrino distributions. There is one file
                       per neutrino flavor, (anti)particle and charged/neutral-current (cc/nc) 
                       interaction permutation.
-- nu*.png           : Example plots showing the data in the corresponding nu*.csv file


Caveats:
--------

The development of the simulations and algorithms for the IceCube Upgrade, as well as the design of 
the detector itself including its geometry, is ongoing. Significant changes are possible, and the 
contents of this data release should be considered preliminary.

A few notable cases are:

- The readout scheme for the mDOMs and DEggs is not yet fully defined. Placeholders for now are: 
    (a) DEgg: similar readout to the IceCube DOM assumed
    (b) 10 ns time resolution assumed (pessimistic)
- Only cascade reconstruction is currently supplied. In the future track+cascade reconstruction 
  ("starting" events) will be utilised, improving the direction reconstruction resolution for events
  with extended muon tracks.
- Background MC (such as atmospheric muons) is not supplied. The veto capabilities of the IceCube 
  Upgrade are expected to significantly exceed IceCube/DeepCore, and the background contamination at 
  analysis level is expected to be low and can be considered neglible to first order.
- Correlated noise across multiple PMTs due to radioactive decays in the glass surrounding the 
  optical modules is not currently modelled.


Conventions:
------------
- The neutrino direction (defined by zenith and azimuth angles) is defined as the direction the 
  neutrino appears to come from (NOT the direction of travel).
- Zenith angle is defined in the interval [0, pi]: 0 = "down-going" (from Southern sky),
  pi = "up-going" (from Northern sky, crossing the Earth). 
- Azimuth angle is defined in the interval [0, 2pi].


events/neutrino_mc.csv:
----------------
This file contains simulated neutrinos with one event per row with the following information:
- true_energy         : Simulated energy [of neutrino, [GeV]
- true_zenith         : Simulated zenith angle of neutrino, [rad]
- true_azimuth        : Simulated azimuth angle of neutrino, [rad]
- reco_energy         : Reconstructed energy of neutrino, [GeV]
- reco_zenith         : Reconstructed zenith angle of neutrino, [rad]
- reco_azimuth        : Reconstructed azimuth angle of neutrino, [rad]
- pdg                 : Code representing particle type, accoridng to the Particle Data Group (PDG)
                        MC particle numbering scheme
- pid                 : Reconstructed particle ID: 0 = cascade, 1 = track
- interaction_type    : Neutrino-ice interaction type: 0 = quasielastic, 1 = resonance, 2 = deep 
                        inelastic scattering, 3 = other (see [2] for details)
- current_type        : Interaction propagator current type: 0 = neutral current, 1 = charged current
- weight              : Weight per event.
                        Expected number of events = flux [m^-2 s^-1] * livetime [s] * weight
- xsec                : Cross section of the neutrino-ice interaction, [1e-38 cm^2]
                        See `EvtXSec` variable in [2] details
- dxsec               : Differential cross section of the neutrino-ice interaction. Units depends on
                        the interaction type.
                        See `EvtDXSec` variable in [2] details, including units.
- x                   : Bjorken x of the neutrino-ice interaction (see [2])
- y                   : Inelasticity of the neutrino-ice interaction (see [2])
- W                   : Invariant hadronic mass of the neutrino-ice interaction, [GeV] (see [2])
- Q2                  : Momentum transfer (Q^2) of the neutrino-ice interaction, [GeV^2] (see [2])

binned/nu*.csv:
----------------
These files contain binned distributions for simulated neutrinos, where the binning is 1D in true 
energy (log10). Each row in the files corresponds to a single bin. The distributions provided are:
- Elow                : Bin lower edge, [GeV].
- Ehigh               : Bin upper edge, [GeV].
- Aeff                : Effective area (all sky), [m^2].
- Ebias               : Median bias of reconstructed neutrino energy w.r.t. true neutrino energy, [GeV].
                        Ebias = median{ Ereco - Etrue }
                        This is an indicator of invisible final state energy (e.g. final state neutrinos).
- |dE|                : Median (absolute) error between reconstructed and true neutrino energy, following 
                        bias subtraction, [GeV].
                        |dE| = median{ |Ereco - Etrue - Ebias| }
                        This is an indicator of spread of reconstructed energies, e.g. resolution.
- dPsi                : Median angular error between reconstructed and true neutrino direction, [deg].
- TrackFraction       : Fraction of events that are identified as track-like (rest are cascade-like).


References:
-----------

[1] A. Ishihara et al. PoS-ICRC2019-1031

[2] C. Andreopoulos et al. arXiv:1510.05494 [hep-ph]

