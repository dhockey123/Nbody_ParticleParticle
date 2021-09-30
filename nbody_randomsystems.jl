using GLMakie

const G = 6.67e-11

Threads.nthreads()
Base.@kwdef mutable struct Params
    dt::Float64
    nbody::Int64
    masses::Array{Float64}

    x_positions::Array{Float64}
    y_positions::Array{Float64}
    z_positions::Array{Float64}
   
    x_velocities::Array{Float64} 
    y_velocities::Array{Float64} 
    z_velocities::Array{Float64} 
   
    x_accelerations::Array{Float64} 
    y_accelerations::Array{Float64} 
    z_accelerations::Array{Float64} 

    Kinetic_Es::Vector{Float64} 
    Potential_Es::Vector{Float64} 
end

dt = 500
nbody = 200
min_lim = -5e3
max_lim = 5e3
axis_lim_min = -1e5
axis_lim_max = 1e5

masses = ones(nbody)
pos_x = rand(min_lim:max_lim, nbody)
pos_y = rand(min_lim:max_lim, nbody)
pos_z = rand(min_lim:max_lim, nbody)
x_velocities = rand(-0.1:0.001:0.1, nbody)
y_velocities = rand(-0.1:0.001:0.1, nbody)
z_velocities = rand(-0.1:0.001:0.1, nbody)
x_accelerations = zeros(Float64, nbody)
y_accelerations = zeros(Float64, nbody)
z_accelerations = zeros(Float64, nbody)

init_params = Params(dt, nbody, masses, 
                    pos_x, pos_y, pos_z, 
                    x_velocities, y_velocities, z_velocities, 
                    x_accelerations, y_accelerations, z_accelerations,
                    zeros(Int64, nbody), zeros(Int64, nbody))

# Gets acceleration for each body by finding the total force projected on the body by all other bodies in the system; O = n²
# softening is applied to prevent large forces as the distance between bodies (r) -> 0. This is not a good idea for specific multi body systems and only worth
# using when bodies are expected to collide
function get_accelerations(par::Params)
    Threads.@threads for i in eachindex(par.x_positions)
        A_x, A_y, A_z = 0,0,0
        softening = 30
        for p in eachindex(par.x_positions)
            if p != i
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 + (par.z_positions[i] - par.z_positions[p])^2+softening^2)
                A_x += -masses[p]*(par.x_positions[i]-par.x_positions[p])/r^3
                A_y += -masses[p]*(par.y_positions[i]-par.y_positions[p])/r^3
                A_z += -masses[p]*(par.z_positions[i]-par.z_positions[p])/r^3
                par.x_accelerations[i] = A_x
                par.y_accelerations[i] = A_y
                par.z_accelerations[i] = A_z
            end
        end
    end
end

# Init velocity for first half time step
function init_velocities(par::Params)
    for i in eachindex(par.x_velocities)
        for p in eachindex(par.x_velocities)
            if p!=i
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 + (par.z_positions[i] - par.x_positions[p])^2)
                Vx = par.x_velocities[i] - par.dt*0.5G*masses[p]*(par.x_positions[i]-par.x_positions[p])/r^1.5
                Vy = par.y_velocities[i] - par.dt*0.5G*masses[p]*(par.y_positions[i]-par.y_positions[p])/r^1.5
                Vz = par.z_velocities[i] - par.dt*0.5G*masses[p]*(par.z_positions[i]-par.z_positions[p])/r^1.5

                par.x_velocities[i] = Vx
                par.y_velocities[i] = Vy
                par.z_velocities[i] = Vz
                
            end
        end
    end
end

# Leapfrog solver for solutions at each time step
function leapfrog(par::Params)
    for i in eachindex(par.x_positions)
        Vx_update = par.x_velocities[i] + par.dt*par.x_accelerations[i]
        Vy_update = par.y_velocities[i] + par.dt*par.y_accelerations[i]
        Vz_update = par.z_velocities[i] + par.dt*par.z_accelerations[i]

        X_pos_update = par.x_positions[i] + par.dt*Vx_update
        Y_pos_update = par.y_positions[i] + par.dt*Vy_update
        Z_pos_update = par.z_positions[i] + par.dt*Vz_update

        par.x_velocities[i] = Vx_update
        par.y_velocities[i] = Vy_update
        par.z_velocities[i] = Vz_update

        par.x_positions[i] = X_pos_update
        par.y_positions[i] = Y_pos_update
        par.z_positions[i] = Z_pos_update
    end
end

fig = Figure(resolution = (1200, 800))
ax1 = Axis3(fig[1,1])

drawx,drawy,drawz = [], [], []
line = []
for i in eachindex(init_params.x_positions)
    append!(drawx, [Node([init_params.x_positions[i]])])
    append!(drawy, [Node([init_params.y_positions[i]])])
    append!(drawz, [Node([init_params.z_positions[i]])])
    append!(line, [Node(Point3f0[( init_params.x_positions[i][], 
                                   init_params.y_positions[i][],
                                   init_params.z_positions[i][] )])])
    lines!(line[i], color=:black)
    scatter!(drawx[i], drawy[i], drawz[i], color = :black,  markersize = init_params.masses[i]*2000)
end

init_velocities(init_params)

function plotme()
    record(fig, "test.gif") do io
        for i in 1:200
            get_accelerations(init_params)
            leapfrog(init_params)

            for i in eachindex(init_params.x_positions)
                drawx[i][] = [init_params.x_positions[i]] 
                drawy[i][] = [init_params.y_positions[i]] 
                drawz[i][] = [init_params.z_positions[i]] 

                # Draws lines to trace body movement.. Comment out from next 6 lines for faster draw times
                new_point = Point3f0(init_params.x_positions[i], init_params.y_positions[i], init_params.z_positions[i])
                if length(line[i][]) > 20
                    line[i][] = push!(line[i][][end-20:end], new_point)
                else
                    line[i][] = push!(line[i][], new_point)
                end

            end
            
            ax1.title = "t = $(round(i*init_params.dt/86400, digits=4)) days"
            recordframe!(io)
        end
    end
end
plotme()

