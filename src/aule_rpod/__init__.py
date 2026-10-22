"""aule_rpod -- docking-window analysis and RPOD approach design.

Layered architecture (each layer depends only on the layers above it):

    config / constants          single source of truth, typed + validated
    timesys / frames            UTC/TT/UT1, IAU 2006/2000A GCRS<->ITRS, LVLH
    ephemeris                   Sun (ERFA epv00; JPL DE421 + analytic for V&V)
    orbit                       SSO design, secular-J2 and Cowell-J2 propagators
    attitude                    STK .a I/O, quaternion conventions, interpolation,
                                closed-form attitude-law identification
    geometry                    ground-station topocentric geometry, eclipse,
                                docking-port illumination
    events                      ms-accurate event detection + interval algebra
    windows                     docking-window logic and near-miss diagnostics
    rpod                        CW guidance, strategy trade, safety, Monte Carlo
    viz / products              figures, Gantt, CSV/JSON products
"""

__version__ = "1.0.0"
