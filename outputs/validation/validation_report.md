| ID | Check | Key result | Verdict |
|---|---|---|---|
| V1 | Attitude law & quaternion convention | fit residual 2.2e-11 rad; SLERP vs law 0.187 deg; NLERP vs law 0.217 deg | PASS |
| V2 | Sun: ERFA vs JPL DE421 / analytic | DE421 0.01" ; almanac 19.0" ; interpolation 4.3e-03" | PASS |
| V3 | GCRS->ITRS vs ERFA direct & astropy | 1.0e-12 rad; ITRS pos 0.020 m; elev 4.9e-07 deg | PASS |
| V4 | Secular J2 vs Cowell J2; LTAN flown | mean reading: +0.0239 deg/day phase drift; osculating reading: +10.98 deg/day; LTAN flown 10.50000-10.50004 h over 285 nodes | PASS |
| V5 | Event-finder grid independence (1/2/5 s) | max diff 6.7e-04 s | PASS |
| V6 | Independent astropy chain, top candidates | longest-span difference <= 0.71 s | PASS |
| V7 | Orekit 13: propagation, elevation, AOS/LOS, eclipse | 0.019 m; 1.3e-06 deg; AOS/LOS 0.91 ms; penumbra 13.8 ms | PASS |
