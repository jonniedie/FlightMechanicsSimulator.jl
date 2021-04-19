"""
    SixDOFAeroEuler
    SixDOFAeroEuler(x::AbstractVector)
    SixDOFAeroEuler(dss::DSState)

Six degrees of freedom dynamic system using TAS, α, β for velocity representation and Euler
angles for attitude.

Flat Earth hypothesis is applied and Earth reference frame is considered inertial.

It is considered that the aircraft xb-zb plane is a plane of symmetry so that Jxy and Jyz
cross-product of inertia are zero and will not be taken into account.

# Fields
- `x::SVector{13, T}`: state vector.
  - [tas (m/s), α (rad), β (rad), ϕ (rad), θ (rad), ψ (rad), p (rad/s), q (rad/s), r (rad/s),
    x (m), y (m), z (m), pow (%)].
"""
struct SixDOFAeroEuler{T}<:DSState{T}
    x::SVector{13, T}
end

SixDOFAeroEuler(x::AbstractVector) = SixDOFAeroEuler(SVector{13, eltype(x)}(x))

SixDOFAeroEuler(dss::DSState) = SixDOFAeroEuler(
    [
        get_tasαβ(dss)...,
        get_euler_angles(dss)[3:-1:1]...,
        get_ang_vel_body(dss)...,
        get_earth_position(dss)...,
        get_engine_power(dss),
    ]
)

get_x_names(dss::SixDOFAeroEuler) = [:tas, :α, :β, :ϕ, :θ, :ψ, :p, :q, :r, :x, :y, :z, :pow]

# Mandatory getters for DSState
get_earth_position(dss::SixDOFAeroEuler) = get_x(dss)[10:12]
get_height(dss::SixDOFAeroEuler) = -get_x(dss)[12]

get_euler_angles(dss::SixDOFAeroEuler) = get_x(dss)[6:-1:4]

get_tas(dss::SixDOFAeroEuler) = get_x(dss)[1]
get_α(dss::SixDOFAeroEuler) = get_x(dss)[2]
get_β(dss::SixDOFAeroEuler) = get_x(dss)[3]
get_tasαβ(dss::SixDOFAeroEuler) = get_x(dss)[1:3]
get_body_velocity(dss::SixDOFAeroEuler) = wind2body(get_tas(dss), 0, 0, get_α(dss), get_β(dss))
get_horizon_velocity(dss::SixDOFAeroEuler) = body2horizon(
    get_body_velocity(dss)..., get_euler_angles(dss)...
)

get_ang_vel_body(dss::SixDOFAeroEuler) = get_x(dss)[7:9]
get_euler_angles_rates(dss::SixDOFAeroEuler) = pqr_2_ψθϕ_dot(
    get_ang_vel_body(dss)..., get_euler_angles(dss)[2:3]...
)

get_engine_power(dss::SixDOFAeroEuler) = get_x(dss)[13]

# Mandatory getters for DSStateDot
get_tas_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[1]
get_α_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[2]
get_β_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[3]
get_tasαβ_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[1:3]
get_accel_body(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = tasαβ_dot_to_uvw_dot(
    get_tasαβ(dssd)..., get_tasαβ_dot(dssd)...
)

get_ang_accel_body(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[7:9]
get_engine_power_dot(dssd::DSStateDot{S, N, T}) where {S<:SixDOFAeroEuler, N, T} = get_xdot(dssd)[13]



function state_eqs(dss::SixDOFAeroEuler, time, mass, inertia, forces, moments, h, pow_dot)

    xdot = sixdof_aero_earth_euler_fixed_mass(
         time,
         get_x(dss),
         mass,
         inertia,
         forces,
         moments,
         h
    )

    xdot = [xdot..., pow_dot]

    return DSStateDot(dss, xdot)
end


"""
    sixdof_aero_earth_euler_fixed_mass!(x_dot, time, x, mass, inertia, forces, moments, h)

Mutating version of [six_dof_aero_euler_fixed_mass](@ref).
"""
function sixdof_aero_earth_euler_fixed_mass!(
    x_dot, time, x, mass, inertia, forces, moments, h,
    )
    # C     x(1)  -> tas (m/s)
    # C     x(2)  -> α (rad)
    # C     x(3)  -> β (rad)
    # C     x(4)  -> ϕ (rad)
    # C     x(5)  -> θ (rad)
    # C     x(6)  -> ψ (rad)
    # C     x(7)  -> p (rad/s)
    # C     x(8)  -> q (rad/s)
    # C     x(9)  -> r (rad/s)
    # C     x(10) -> North (m)
    # C     x(11) -> East (m)
    # C     x(12) -> Altitude (m)

    # Assign state
    @log tas = x[1]
    @log α = x[2]
    @log β = x[3]
    @log ϕ = x[4]
    @log θ = x[5]
    @log ψ = x[6]
    @log p = x[7]
    @log q = x[8]
    @log r = x[9]

    # Unpack forces
    @log Fx, Fy, Fz = forces
    # Unpack moments
    @log L, M, N = moments
    # Unpack angular momentum contributions
    @log hx, hy, hz = h
    # Unpack inertia
    Ixx, Iyy, Izz = inertia[1, 1], inertia[2, 2], inertia[3, 3]
    Ixz = inertia[1, 3]

    # Get ready for state equations
    @log u, v, w = wind2body(tas, 0, 0, α, β)

    sψ, cψ = sin(ψ), cos(ψ)
    sθ, cθ = sin(θ), cos(θ)
    sϕ, cϕ = sin(ϕ), cos(ϕ)

    # Force equations
    udot = r * v - q * w + Fx / mass
    vdot = p * w - r * u + Fy / mass
    wdot = q * u - p * v + Fz / mass

    x_dot[1], x_dot[2], x_dot[3] = uvw_dot_to_tasαβ_dot(u, v, w, udot, vdot, wdot)

    # Kinematics
    x_dot[6], x_dot[5], x_dot[4] = pqr_2_ψθϕ_dot(p, q, r, θ, ϕ)

    # Moments
    pq = p * q
    qr = q * r

    rhy_qhz = (r * hy - q * hz)
    qhx_phy = (q * hx - p * hy)

    # If inertia is constant this terms are constant too.
    # TODO: think about passing them as arguments to improve speed.
    IxzS = Ixz^2
    xpq = Ixz * (Ixx - Iyy + Izz)
    gam = Ixx * Izz - IxzS
    xqr = Izz * (Izz - Iyy) + IxzS
    zpq = (Ixx - Iyy) * Ixx + IxzS
    ypr = Izz - Ixx

    x_dot[7] = (xpq * pq - xqr * qr + Izz * (L + rhy_qhz) + Ixz * (N + qhx_phy)) / gam
    x_dot[8] = (ypr * p * r - Ixz * (p^2 - r^2) + M - r * hx + p * hz) / Iyy
    x_dot[9] = (zpq * pq - xpq * qr + Ixz * (L + rhy_qhz) + Ixx * (N + qhx_phy)) / gam

    # Navigation
    t1 = sϕ * cψ
    t2 = cϕ * sθ
    t3 = sϕ * sψ
    s1 = cθ * cψ
    s2 = cθ * sψ
    s3 = t1 * sθ - cϕ * sψ
    s4 = t3 * sθ + cϕ * cψ
    s5 = sϕ * cθ
    s6 = t2 * cψ + t3
    s7 = t2 * sψ - t1
    s8 = cϕ * cθ

    x_dot[10] = u * s1 + v * s3 + w * s6  # North speed
    x_dot[11] = u * s2 + v * s4 + w * s7  # East speed
    x_dot[12] = -u * sθ + v * s5 + w * s8  # Vertical speed (possitive downwards)

end


"""
    sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)

State equations for `SixDOFAeroEuler` dynamic system.

# Arguments

- `time::Number`: time (s).
- `x::AbstractArray`: state vector.
  - [tas (m/s), α (rad), β (rad), ϕ (rad), θ (rad), ψ (rad), p (rad/s), q (rad/s), r (rad/s),
    x (m), y (m), z (m), pow (%)].
- `mass::Number`: mass (kg).
- `inertia::AbstractMatrix`: 3x3 inertia tensor (kg·m^2).
- `forces::AbstractVector`: body axis forces [Fx, Fy, Fz] (N).
- `moments::AbstractVector`: body axis moments (L, M, N) (N·m).
- `h::AbstractVector`:  Additional angular momentum contributions such as those coming from
 spinning rotors (kg·m²/s).

# Returns
- `x_dot::Array`: time derivative of the state vector.

# Notes

It is considered that the aircraft xb-zb plane is a plane of symmetry so that Jxy and Jyz
cross-product of inertia are zero and will not be taken into account.

# See also

[six_dof_aero_euler_fixed_mass!](@ref)
"""
function sixdof_aero_earth_euler_fixed_mass(time, x, mass, inertia, forces, moments, h)
    x_dot = Array{eltype(x)}(undef, 12)
    sixdof_aero_earth_euler_fixed_mass!(x_dot, time, x, mass, inertia, forces, moments, h)
    return x_dot
end
