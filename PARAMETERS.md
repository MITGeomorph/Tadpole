# Tadpole 2.0 parameter reference

All parameters are fields of the struct `p` passed to `Tadpole(init,p)`. "Required" parameters must be set in the run script; the others have the default shown. Grids ("scalar or grid") can be a scalar or a matrix of size `[p.Ny p.Nx]`. Units: m, yr.

## Grid and time

| Parameter | Default | Description |
|---|---|---|
| `p.Nx`, `p.Ny` | required | Number of grid points in x (columns) and y (rows). `init` must be `p.Ny` x `p.Nx`. |
| `p.dx`, `p.dy` | required | Grid spacing (m). |
| `p.dtmax` | required | Time step (yr). The implicit solvers are stable for any time step, so choose it for accuracy. With the explicit fluvial solver, the maximum time step unless `p.fixedTimeStep = 1`. |
| `p.tf` | required | Duration of the run (yr). |
| `p.fixedTimeStep` | 0 | 1 = every time step is `p.dtmax`. 0 = the explicit fluvial solver may shorten a time step to keep its Courant number below `p.StabilityFactor`. Matters only with `p.fluvialSolver = 'explicit'`; the implicit solvers always take steps of `p.dtmax`. With 1 and the explicit solver, Tadpole warns if the Courant number exceeds 1. |
| `p.StabilityFactor` | 0.5 | Advanced. Maximum Courant number of the explicit fluvial solver when `p.fixedTimeStep = 0`. |

## Boundaries, sources and flow routing

| Parameter | Default | Description |
|---|---|---|
| `p.bdy.left`, `.right`, `.upper`, `.lower` | required | Boundary condition on each edge: `'fixed'` (constant elevation; water, sediment and soil leave the model), `'mirror'` (reflection: nothing crosses the edge, a drainage divide), or `'periodic'` (joined to the opposite edge; set both edges of a pair). |
| `p.F` | none | Grid of fixed points: cells with `p.F == 1` keep their elevation and act like fixed boundaries (for example the sea around an island). |
| `p.U` | required | Rock uplift rate relative to the fixed boundaries (m/yr), scalar or grid. |
| `p.P` | required | Mean annual precipitation rate (m/yr), scalar or grid. Water discharge is P averaged over the upstream area times drainage area. |
| `p.routing` | required | Flow routing: `'D8'` (steepest descent), `'Dinf'` (D-infinity, Tarboton 1997) or `'Dms'` (multiple flow directions). |
| `p.mfdexp` | 1 | `'Dms'` only: flow is divided among downslope neighbours in proportion to slope^`p.mfdexp`. Larger values concentrate flow in the steepest direction. |
| `p.flood` | required | 1 = fill closed depressions (Priority-Flood) so that flow crosses lakes to their outlets; 0 = flow ends in depressions. Must be 1 with sediment. |

## River incision

| Parameter | Default | Description |
|---|---|---|
| `p.doStreamPower` | required | 1 = river incision on, 0 = off. |
| `p.Kr` | required | Bedrock erodibility in E = Kr (P A)^m S^n, units m^(1-3m) yr^(m-1). Scalar or grid. |
| `p.m` | required | Discharge exponent m. |
| `p.n` | required | Slope exponent n. |
| `p.Ecr` | required with `'smooth'` | Threshold erosion rate for bedrock (m/yr) of the smoothed threshold F(x) = x - Ecr (1 - exp(-x/Ecr)), x = Kr (P A)^m S^n. Cells where Kr (P A)^m S^n > Ecr count as channels. 0 = no threshold. Rule of thumb: for channel heads at drainage area Ah and slope Sh, Ecr ~ Kr (P Ah)^m Sh^n. |
| `p.channelThreshold` | `'smooth'` | `'smooth'`: every cell follows the law with the smoothed thresholds `p.Ecr`, `p.Ecs`. `'area'`: rivers erode only where drainage area > `p.Ac` (the rule of Tadpole 1); much faster with sediment, but hollows near channel heads can fill with soil. A run script that sets `p.Ac` but not `p.Ecr` runs with `'area'` and a warning. |
| `p.Ac` | required with `'area'` | Drainage area threshold (m^2). |
| `p.fluvialSolver` | `'implicit'` | `'implicit'` (backward time, stable for any time step, required for sediment) or `'explicit'` (forward time, time step limited by the Courant number). |
| `p.pitLimiter` | 1 | Implicit solver: 1 = a cell cannot be lowered below its lowest receiver in one step. |
| `p.doChannelDiffusion` | required | 1 = hillslope transport also acts in channel cells; 0 = channel cells are sinks for hillslope soil (with sediment, that soil joins the river's sediment load). |

## Sediment

Used only with `p.doSediment = 1`, which requires `p.flood = 1` and `p.fluvialSolver = 'implicit'`.

| Parameter | Default | Description |
|---|---|---|
| `p.doSediment` | required | 1 = sediment entrainment and deposition (SPACE; Shobe et al., 2017) on, 0 = off. |
| `p.Ks` | required | Sediment erodibility (units of `p.Kr`). Scalar or grid. |
| `p.Ecs` | required with `'smooth'` | Threshold erosion rate for sediment entrainment (m/yr), as `p.Ecr` but with `p.Ks`. |
| `p.Hstar` | required | Bed roughness length scale H* (m): bedrock is partly exposed where the sediment is thinner than about H*. |
| `p.Vs` | required | Effective settling velocity (m/yr), as in Davy & Lague (2009). Scalar or grid. Larger values (relative to `p.P`) deposit more sediment and need more solver iterations. |
| `p.phi` | 0 | Porosity of sediment and soil. |
| `p.Ff` | 0 | Fraction of eroded bedrock that becomes wash load and leaves the model. |
| `p.H0` | 0 | Initial sediment thickness (m), scalar or grid. |
| `p.sedTol` | 1e-3 max(\|U\|) dtmax | Convergence tolerance of the sediment solver (m). Must be set if there is no uplift anywhere: about 1e-3 times the typical elevation change per time step. |
| `p.sedMaxIter` | 200 | Advanced. Maximum number of solver iterations per time step. |

## Hillslopes

| Parameter | Default | Description |
|---|---|---|
| `p.doNLDiffusion` | required | 1 = hillslope transport on, 0 = off. |
| `p.Knl` | required | Hillslope transport coefficient (m^2/yr). Scalar or grid. |
| `p.Sc` | required | Critical slope of the nonlinear transport law. `Inf` gives linear diffusion. |
| `p.hillSolver` | `'auto'` | `'auto'`: for each time step the cheaper of `'explicit'` (forward time, divided into stable sub-steps) and `'implicit'` (backward Euler, stable for any time step) for the same law; with `p.Sc = Inf`, an implicit linear solver. `'explicit'` or `'implicit'` forces one. |
| `p.hillPhiSwitch` | 0.3 | Advanced. `'auto'`: the explicit solver is used when dt Knl/dx^2 <= `p.hillPhiSwitch`, the implicit solver above. |
| `p.hillStabilityFactor` | 0.8 | Advanced. The explicit hillslope solver takes sub-steps of this fraction of its maximum stable time step, 1/(2 D (1/dx^2 + 1/dy^2)), with D the largest diffusivity on the grid. Must be in (0, 1]. |
| `p.hillN` | 5 | Advanced. Number of terms after the first in the series approximation of the nonlinear diffusivity (Ganti et al., 2012). |

## Plotting and output

| Parameter | Default | Description |
|---|---|---|
| `p.doDrawPlot` | required | 1 = display the surface during the run. |
| `p.plotint` | required with plotting | Redraw the plot every `p.plotint` time steps. |
| `p.plottype` | required with plotting | `'elevation'`, `'hillshade'`, `'sediment'` (shaded relief colored by sediment thickness; requires sediment) or `'perspective'` (oblique 3-D view). |
| `p.plotVE` | automatic | `'perspective'`: vertical exaggeration. |
| `p.doSaveOutput` | required | 1 = save the run to `<p.runname>.mat` and an image of the final surface to `<p.runname>.png`. |
| `p.saveint` | required with saving | Save the elevation grid every `p.saveint` time steps. |
| `p.runname` | `'default_run_name'` | Name of the output files. |

## Outputs

`solution = Tadpole(init,p)` returns the final elevation grid, or, with `p.doSaveOutput = 1`, all saved grids (`p.Ny` x `p.Nx` x number of saves).

With `p.doSaveOutput = 1`, `<p.runname>.mat` contains:
- `elevation`: elevations (`p.Ny` x `p.Nx` x number of saves), including the initial and final surfaces;
- `t`: the time of each saved grid (yr);
- `p`: the parameters;
- `sediment` (with sediment): sediment thickness at the same times; bedrock elevation is `elevation - sediment`.

`TadpoleMovie(elevation,p,t,frame_interval,basename[,sediment])` writes one .png image for every `frame_interval`-th saved grid, using `p.plottype`, into a new subfolder named `basename`.
