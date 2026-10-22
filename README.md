# Tadpole 2.0

MIT landscape evolution model

Tadpole simulates the evolution of topography on a regular grid by rock uplift, river incision, sediment transport and deposition by rivers, and hillslope soil transport. It runs in MATLAB, with the computationally demanding parts written as C MEX files. Version 2.0 is a major revision of the 2017 (and earlier) version of Tadpole: the numerical methods for flow routing, river incision and hillslope transport are new, and rivers can now carry and deposit sediment. See [What's new in Tadpole 2.0](#whats-new-in-tadpole-20).

A parameter reference is in [PARAMETERS.md](PARAMETERS.md).


## Getting started

1. **Download** Tadpole: clone this repository, or use *Code > Download ZIP* on GitHub.
2. **Compile the MEX files.** Tadpole needs a C compiler that MATLAB supports: on Windows, Microsoft Visual C++ (free with Visual Studio Community) or the MinGW-w64 add-on from MATLAB's Add-On Explorer; on macOS, Xcode or its command line tools (`xcode-select --install`); on Linux, gcc. In MATLAB, `mex -setup C` shows or chooses the compiler. Then, from the Tadpole folder:
   ```matlab
   compile_mex
   ```
   This builds six MEX files for your platform (`.mexw64`, `.mexmaca64`, `.mexmaci64` or `.mexa64`). You need to do it once, and again only if you change a `.c` or `.h` file.
3. **Run an example.** The example scripts are in the Tadpole folder. Run one from there, or copy it and edit the copy:
   ```matlab
   EXAMPLE_RUNSCRIPT_valleys
   ```

| Example | What it shows | Run time* |
|---|---|---|
| `EXAMPLE_RUNSCRIPT_island.m` | Erosion of a volcanic island with no uplift. Irregular boundaries: fixed points inside the domain (the sea). Precipitation that varies in space, a smoothed erosion threshold, and channels as sinks for hillslope soil. | 2 min |
| `EXAMPLE_RUNSCRIPT_rivers.m` | River networks cutting into an uplifting plateau on a large grid (1000 x 1000), with fixed boundaries all around and river incision only. Long time steps (100,000 years) with the implicit solver. | 40 s |
| `EXAMPLE_RUNSCRIPT_valleys.m` | Evenly spaced ridges and valleys on a grid of 10 m cells: nonlinear hillslope transport, including in channels, and river incision without thresholds. Periodic and fixed boundaries; saving output and making movie frames. | 2 min |
| `EXAMPLE_RUNSCRIPT_sediment.m` | A mountain range shedding sediment onto a foreland, after Yuan et al. (2019): alluvial fans that merge into a plain crossed by shifting channels. Sediment entrainment and deposition, uplift that varies in space, multiple-flow routing, mirror and periodic boundaries. | 3.5 min |

\*On a Windows laptop without plotting; plotting adds time.

Every run script has the same structure: it sets the parameters in a struct `p`, makes an initial elevation grid `init`, and calls `solution = Tadpole(init,p)`. Each parameter is described in the script where it is set and in [PARAMETERS.md](PARAMETERS.md).


## What the model solves

Elevation z(x, y, t) changes by rock uplift, river incision and hillslope transport,

    dz/dt = U - E + (hillslope transport)

- **Rock uplift** at rate U relative to the fixed boundaries (`p.U`, scalar or grid).
- **River incision** by the stream power law, E = F(Kr Ω), with Ω = (P A)^m S^n, where P A is the water discharge (precipitation rate P times drainage area A, with P averaged over the upstream area), S the slope along the flow paths, and F a smoothed erosion threshold (Shobe et al., 2017):

      F(x; Ec) = x - Ec (1 - exp(-x/Ec))

  F(x; 0) = x (no threshold). F is close to x²/(2 Ec) for x << Ec and to x - Ec for x >> Ec, so erosion fades out gradually below the threshold rate Ec instead of switching off. Alternatively, rivers erode only where A exceeds a threshold area Ac (`p.channelThreshold = 'area'`).
- **Sediment** (optional, `p.doSediment = 1`): the SPACE model of Shobe et al. (2017). Rivers erode bedrock at rate Er = F(Kr Ω; Ecr) exp(-H/H*), entrain sediment at rate Es = F(Ks Ω; Ecs) (1 - exp(-H/H*)), and deposit sediment at rate D = Vs Qs/Q, where H is the sediment thickness, H* a bed roughness length scale, Vs an effective settling velocity, Qs the sediment flux and Q the water discharge. Lakes trap the sediment that reaches them until they fill to their spill level. Hillslope soil delivered to channels joins the rivers' sediment load.
- **Hillslope transport** by the nonlinear law of Roering et al. (1999), q = -K ∇z / (1 - (|∇z|/Sc)²), with the diffusivity approximated by the series of Ganti et al. (2012), or linear diffusion (`p.Sc = Inf`). With sediment, hillslopes are soil-mantled: soil is produced from bedrock where hillslopes erode.

Each time step applies uplift, then hillslope transport, then river incision and sediment transport (operator splitting).


## Numerical methods

- **Flow routing** by steepest descent (D8), D-infinity (Tarboton, 1997) or multiple flow directions in proportion to slope raised to a power (Freeman, 1991). Closed depressions are filled by Priority-Flood+ε (Barnes et al., 2014) on every time step, so flow crosses lakes to their outlets; the flow weights are computed on the filled surface, so every flow path leads downhill to a boundary or fixed point. Cells are ordered with each cell before its donors by a topological sort, which works for single and multiple flow directions.
- **River incision and sediment** are solved implicitly along that order (Braun & Willett, 2013), so the time step is not limited by a Courant number. With sediment, the coupling between deposition and the sediment flux from upstream is solved by Gauss-Seidel iteration (Yuan et al., 2019), with an exact mass budget. An explicit solver (`p.fluvialSolver = 'explicit'`) remains available without sediment.
- **Hillslope transport** uses the conservative finite-difference scheme of Andrews & Hanks (1985) as described by Perron (2011). By default (`p.hillSolver = 'auto'`), each time step uses the cheaper of two solvers for the same law: explicit sub-steps when the time step is short relative to dx²/K, and a backward-Euler step solved by Picard iteration with Anderson acceleration (Walker & Ni, 2011) when it is long. Linear diffusion is solved implicitly.

The model works for any combination of fixed, mirror and periodic boundaries, and for fixed points inside the domain (`p.F`).


## What's new in Tadpole 2.0

The previous version on GitHub dates from 2017. Since then Tadpole has changed in two stages: unpublished development by Taylor Perron (2017-2026), and a 2026 overhaul that included a rewrite of the flow routing and depression handling and the addition of sediment and new solvers.

### Flow routing and depressions

The flooding code of the previous version had several problems, which this revision fixes.

- **Old method.** Depressions were handled by finding each one's lowest pour point and merging basins, with a fixed limit on the number of neighbouring basins (1000; rough surfaces could stop the run with "Memory overrun"). Flow directions for cells outside depressions came from the unfilled surface and those inside from the filled one. At a depression's spill point that combination could send flow into the depression and back out again (a flow loop), which made drainage areas wrong in rare cases. Drainage area was accumulated recursively, which could overflow the C stack on large grids with long flow paths.
- **New method.** Depressions are filled by Priority-Flood+ε (Barnes et al., 2014), for all combinations of fixed, mirror and periodic boundaries, with a minimum increment that keeps the filled slopes numerically robust. Flow weights for all cells are computed on the filled surface, so every transfer of flow goes strictly downhill and no loops are possible. Drainage area and discharge are accumulated without recursion, along a topological ordering of the cells (Kahn's algorithm; Braun and Willett, 2013), which works for D8, D-infinity and multiple flow directions.
- **New outputs** of the routing functions, used by the new solvers: the cell ordering (`g.Stack`), the receivers of each cell (`g.Receiver`), the filled surface (`g.Zf`) and a label for each lake (`g.Lake`). All three routing functions (`mexD8`, `mexDinf`, `mexDms`) share one interface.
- **Other fixes:** drainage area at corners between two mirror boundaries (0.1-0.5% too large before); slopes of cells on mirror boundaries; a neighbour that the multiple-flow routing skipped. 'Dms' routing has a slope exponent (`p.mfdexp`).
- Rivers that flow into a lake erode down to the lake surface, not to the lake floor, because stream power depends on the slope of the water surface.

### River incision

- **Implicit solver** (`p.fluvialSolver = 'implicit'`, the default), in C: stable for any time step.
- **Erosion law** E = Kr (P A)^m S^n, with water discharge P A instead of drainage area, so precipitation can vary in space (`p.P` as a grid). 
- **Smoothed erosion thresholds** (`p.Ecr`, `p.Ecs`, in m/yr) instead of an area threshold. With an area threshold (Ac), a hollow near a channel head whose drainage area drops below Ac stops eroding, and hillslope soil fills it; a threshold on stream power responds to the fill, because filling steepens the hollow. The area threshold is still available (`p.channelThreshold = 'area'`, `p.Ac`) and is much faster with sediment.

### Sediment (new)

- The SPACE model (Shobe et al., 2017): bedrock and sediment layers, entrainment and deposition, porosity (`p.phi`) and wash load (`p.Ff`). Solved implicitly with Gauss-Seidel iteration (Yuan et al., 2019), with an exact mass budget each time step (`g.Budget`).
- Lakes trap sediment and fill up to their spill level.
- Hillslope soil delivered to channels enters the rivers' sediment load.
- Outputs: sediment thickness `g.H`, bedrock elevation `g.R`, sediment flux `g.Qs`; saved runs contain the sediment thickness (`sediment`); plot type 'sediment'.

### Hillslopes

- Nonlinear transport with a critical slope (Roering et al., 1999) replaces linear diffusion by ADI and optional landslides. Linear diffusion is the case `p.Sc = Inf`.
- Explicit, sub-stepped solver and an implicit Picard solver in C for the same law, chosen automatically for each time step (`p.hillSolver = 'auto'`).
- Channel cells can be sinks for hillslope soil (`p.doChannelDiffusion = 0`). Faces between a diffusing cell and a sink use the full diffusivity of the diffusing cell (the flux into sinks was halved before), and a channel cell passes no hillslope flux to a lower hillslope neighbour.

### Other changes

- Plot types 'elevation', 'hillshade', 'sediment' and 'perspective' (with `p.plotVE`); `TadpoleMovie` makes movie frames from saved runs.
- `compile_mex.m` builds the MEX files on Windows, macOS (Intel and Apple silicon) and Linux. Compiled MEX files are no longer distributed.
- Removed: landslides, linear diffusion by ADI, the plot types 'mesh', 'contour', 'drainage area', 'curvature', 'shade' and 'color shade', and the old examples.


## Moving a run script from an older (pre-2026) version of Tadpole

Start from one of the new example scripts rather than editing an old script. The main differences:

| Legacy Tadpole | Tadpole 2.0 |
|---|---|
| `p.E` (uplift rate) | `p.U` |
| `p.w` (slope exponent) | `p.n` (converted automatically, with a warning) |
| `p.Kf`, `p.Kw`, `p.wexp` (E = Kf Kw A^wexp/dx A^m S^n) | `p.Kr` (E = Kr (P A)^m S^n); different units and values |
| `p.thetac` (threshold on A^m S^n) | `p.Ecr` (threshold erosion rate, m/yr), or `p.channelThreshold = 'area'` with `p.Ac` |
| `p.doDiffusion`, `p.D` (linear diffusion) | `p.doNLDiffusion`, `p.Knl`, with `p.Sc = Inf` for linear diffusion |
| `p.doLandslides`, `p.Sc` | nonlinear transport with critical slope `p.Sc` |
| `p.Courant` | `p.StabilityFactor` (explicit solvers only) |
| `p.doAdaptiveTimeStep` | `p.fixedTimeStep`, with the opposite meaning (converted automatically, with a warning) |
| `p.zfactor`, plot types | `p.plotVE`; plot types 'elevation', 'hillshade', 'sediment', 'perspective' |
| no precipitation | `p.P` (m/yr), scalar or grid |

`Tadpole` stops with an error that points here if a run script sets `p.doDiffusion`, `p.D`, `p.doLandslides`, `p.thetac` or `p.Courant`.


## Citing Tadpole

If you use Tadpole for work that results in a publication, please acknowledge Tadpole in the text and cite one or more of the following papers, as appropriate:

Perron, J.T., W.E. Dietrich and J.W. Kirchner (2008), Controls on the spacing of first-order valleys. Journal of Geophysical Research, 113, F04016, doi: 10.1029/2007JF000977.

Perron, J.T., J.W. Kirchner and W.E. Dietrich (2009), Formation of evenly spaced ridges and valleys. Nature, 460, 502–505, doi: 10.1038/nature08174.

Perron, J.T., P.W. Richardson, K.L. Ferrier, and M. Lapôtre (2012), The root of branching river networks. Nature, 492, 100–103, doi: 10.1038/nature11672.

Richardson, P.W., J.T. Perron, S.R. Miller and J.W. Kirchner (2020), Modeling the formation of topographic asymmetry by aspect-dependent erosional processes and lateral channel migration. Journal of Geophysical Research: Earth Surface, 125, e2019JF005377, doi: 10.1029/2019JF005377.

Richardson, P.W., J.T. Perron, S.R. Miller and J.W. Kirchner (2020), Unraveling the mysteries of asymmetric topography at Gabilan Mesa, California. Journal of Geophysical Research: Earth Surface, 125, e2019JF005378, doi: 10.1029/2019JF005378.

Please also cite the papers describing the methods you use (see the references below), in particular Shobe et al. (2017) if you use sediment.


## References

Andrews, D.J., and Hanks, T.C. (1985), Scarp degraded by linear diffusion: Inverse solution for age. Journal of Geophysical Research, 90, 10193–10208.

Barnes, R., Lehman, C., and Mulla, D. (2014), Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences, 62, 117–127.

Braun, J., and Willett, S.D. (2013), A very efficient O(n), implicit and parallel method to solve the stream power equation governing fluvial incision and landscape evolution. Geomorphology, 180–181, 170–179.

Davy, P., and Lague, D. (2009), Fluvial erosion/transport equation of landscape evolution models revisited. Journal of Geophysical Research, 114, F03007.

Freeman, T.G. (1991), Calculating catchment area with divergent flow based on a regular grid. Computers & Geosciences, 17(3), 413–422.

Ganti, V., Passalacqua, P., and Foufoula-Georgiou, E. (2012), A sub-grid scale closure for nonlinear hillslope sediment transport models. Journal of Geophysical Research, 117, F02012.

Perron, J.T. (2011), Numerical methods for nonlinear hillslope transport laws. Journal of Geophysical Research, 116, F02021.

Roering, J.J., Kirchner, J.W., and Dietrich, W.E. (1999), Evidence for nonlinear, diffusive sediment transport on hillslopes and implications for landscape morphology. Water Resources Research, 35, 853–870.

Shobe, C.M., Tucker, G.E., and Barnhart, K.R. (2017), The SPACE 1.0 model: a Landlab component for 2-D calculation of sediment transport, bedrock erosion, and landscape evolution. Geoscientific Model Development, 10, 4577–4604.

Tarboton, D.G. (1997), A new method for the determination of flow directions and upslope areas in grid digital elevation models. Water Resources Research, 33, 309–319.

Walker, H.F., and Ni, P. (2011), Anderson acceleration for fixed-point iterations. SIAM Journal on Numerical Analysis, 49, 1715–1735.

Yuan, X.P., Braun, J., Guerit, L., Rouby, D., and Cordonnier, G. (2019), A new efficient method to solve the stream power law model taking into account sediment deposition. Journal of Geophysical Research: Earth Surface, 124, 1346–1365.


## License

MIT License; see [LICENSE](LICENSE).
