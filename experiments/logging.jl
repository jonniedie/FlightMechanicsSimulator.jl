using FlightMechanicsSimulator
using FlightMechanicsUtils
using OrdinaryDiffEq
using SimulationLogs
using Test

x_stev = [
    502.0 * FT2M,
    0.2392628,
    0.0005061803,
    1.366289,
    0.05000808,
    0.2340769,
    -0.01499617,
    0.2933811,
    0.06084932,
    0.0 * FT2M,
    0.0 * FT2M,
    0.0 * FT2M,
    64.12363,
]
controls_stev = [0.8349601, -1.481766, 0.09553108, -0.4118124]
xcg = 0.35
dssd, outputs = f(
    time,
    SixDOFAeroEuler(x_stev),
    controls_stev,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere(x_stev[12]),
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
)

x_dot = get_xdot(dssd)

# RETRIM to refine flying condition
dssd, controls_trim, outputs_trim, cost = trim(
    SixDOFAeroEuler(x_stev),
    controls_stev,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere(x_stev[12]),
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M),
    0.0,
    0.3,
)

x_trim = get_x(dssd)
x_dot_trim = get_xdot(dssd)

dt = 0.01  # s
t0 = 0.0  # s
t1 = 180.0  # s

x = x_trim
# Transform to Input for simulate
controls = ConstantInput.(controls_trim)

results = FlightMechanicsSimulator._simulate(
    t0,
    t1,
    SixDOFAeroEuler(x),
    controls,
    F16(F16Stevens.MASS, F16Stevens.INERTIA, xcg),
    F16StevensAtmosphere,
    LHDownGravity(FlightMechanicsSimulator.F16Stevens.GD*FT2M);
    solver=RK4(),
    solve_args=Dict(:reltol=>1e-10, :saveat=>dt),
    )

out = get_log(results)

@testset "Logging" begin
    for (i, name) in enumerate(get_x_names(dssd)[1:9])
        @test results[i,:] == getproperty(out, name)
    end
end

out.Fx
