using GLMakie

const G = 6.67e-11

min_lim = -5e5
max_lim = 5e5
axis_lim_min = -2e5
axis_lim_max = 2e5

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

## Two body system

nbody = 2
dt = .05
masses = [3e15, 3e15]
pos_x = [-110, 110]
pos_y = [0,0]
pos_z = [0,0]
x_velocities = [0,0]
y_velocities = [12.65, -12.65]
z_velocities = [-13, +13]
x_accelerations = zeros(Float64, 2)
y_accelerations = zeros(Float64, 2)
z_accelerations = zeros(Float64, 2)

## Five body system
# Uncomment this section for 5-body system, the 2-body will get overwritten

# nbody = 5
# dt = .05
# masses = [3e12, 3e12, 1e11, 1e11, 1e15]
# pos_x = [-110, 110, 0, 0, 0]
# pos_y = [0,0, 50, -50, 0]
# pos_z = [0,0, 0, 0, 0 ]
# x_velocities = [0,0, 46.5, -46.5, 0]
# y_velocities = [12.65, -12.65, 0, 0,0]
# z_velocities = [-13, +13, -12 , +12,0]
# x_accelerations = zeros(Float64, 5)
# y_accelerations = zeros(Float64, 5)
# z_accelerations = zeros(Float64, 5)


init_params = Params(dt, nbody, masses, 
                    pos_x, pos_y, pos_z, 
                    x_velocities, y_velocities, z_velocities, 
                    x_accelerations, y_accelerations, z_accelerations,
                    zeros(Int64, nbody), zeros(Int64, nbody))

###################################################################################################################################
# Gets acceleration for each body by finding the total force projected on the body by all other bodies in the system; O = n²
# softening is applied to prevent large forces as the distance between bodies (r) -> 0. This is not a good idea for specific multi body systems and only worth
# using when bodies are expected to collide
function get_accelerations(par::Params)
    Threads.@threads for i in eachindex(par.x_positions)
        A_x, A_y, A_z = 0,0,0
        par.Potential_Es[i] = 0
        softening = 0
        par.Kinetic_Es[i] = 0
        for p in eachindex(par.x_positions)
            if p != i
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 + (par.z_positions[i] - par.z_positions[p])^2 + softening^2)
                A_x += -par.masses[p]*G*(par.x_positions[i]-par.x_positions[p])/r^3
                A_y += -par.masses[p]*G*(par.y_positions[i]-par.y_positions[p])/r^3
                A_z += -par.masses[p]*G*(par.z_positions[i]-par.z_positions[p])/r^3
                par.Potential_Es[i] += -G*masses[p]*masses[i] / r
                par.x_accelerations[i] = A_x
                par.y_accelerations[i] = A_y
                par.z_accelerations[i] = A_z
            end
        end
    end
end
# Init velocity for first half time step
function init_velocities(par::Params)
    dt = par.dt
    for i in eachindex(par.x_velocities)
        for p in eachindex(par.x_velocities)
            if p!=i
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 + (par.z_positions[i] - par.x_positions[p])^2)
                Vx = par.x_velocities[i] - dt*0.5G*masses[p]*(par.x_positions[i]-par.x_positions[p])/r^3
                Vy = par.y_velocities[i] - dt*0.5G*masses[p]*(par.y_positions[i]-par.y_positions[p])/r^3
                Vz = par.z_velocities[i] - dt*0.5G*masses[p]*(par.z_positions[i]-par.z_positions[p])/r^3

                par.x_velocities[i] = Vx
                par.y_velocities[i] = Vy
                par.z_velocities[i] = Vz
                                 
            end
        end
    end
end

# Leapfrog solver for solutions at each time step
function leapfrog(par::Params)
    dt = par.dt
    for i in eachindex(par.x_positions)
        Vx_update = par.x_velocities[i] + dt*par.x_accelerations[i]
        Vy_update = par.y_velocities[i] + dt*par.y_accelerations[i]
        Vz_update = par.z_velocities[i] + dt*par.z_accelerations[i]

        X_pos_update = par.x_positions[i] + dt*Vx_update
        Y_pos_update = par.y_positions[i] + dt*Vy_update
        Z_pos_update = par.z_positions[i] + dt*Vz_update  

        par.x_velocities[i] = Vx_update
        par.y_velocities[i] = Vy_update
        par.z_velocities[i] = Vz_update

        par.Kinetic_Es[i] += par.masses[i]*(Vx_update^2+Vy_update^2+Vz_update^2)

        par.x_positions[i] = X_pos_update
        par.y_positions[i] = Y_pos_update
        par.z_positions[i] = Z_pos_update
    end
end

fig = Figure(resolution = (1200, 800))
ax1 = Axis3(fig[1:3,1:2])
ax2 = Axis(fig[1,3], title="Energy")
ax3 = Axis(fig[2,3], title="Kinetic E")
ax4 = Axis(fig[3,3], title="Potential E")

drawx,drawy,drawz = [], [], []
line = []
for i in eachindex(init_params.x_positions)
    append!(drawx, [Node([init_params.x_positions[i]])])
    append!(drawy, [Node([init_params.y_positions[i]])])
    append!(drawz, [Node([init_params.z_positions[i]])])
    append!(line, [Node(Point3f0[( init_params.x_positions[i][], 
                                   init_params.y_positions[i][],
                                   init_params.z_positions[i][] )])])
    lines!(ax1, line[i], color=:blue)
    scatter!(ax1, drawx[i], drawy[i], drawz[i], color = :black,  markersize = 3000)
end
xlims!(ax1, -200, 200)
ylims!(ax1, -200, 200)
zlims!(ax1, -200, 200)
ke = Node(Point2f0[(0,init_params.Kinetic_Es[1])])
pe = Node(Point2f0[(0,init_params.Potential_Es[1])])
energy = Node(Point2f0[(0,init_params.Potential_Es[1]+init_params.Kinetic_Es[1])])

lines!(ax2, energy, color=:black)   
lines!(ax3, ke)
lines!(ax4, pe)

init_velocities(init_params)

function plotme()
    record(fig, "test.mp4") do io
        for i in 1:100000
            
            get_accelerations(init_params)
            leapfrog(init_params)
            
            ke_up = Point2f0(i, init_params.Kinetic_Es[1])
            pe_up = Point2f0(i, init_params.Potential_Es[1])
            e_up = Point2f0(i, init_params.Potential_Es[1]+init_params.Kinetic_Es[1])
            ke[] = push!(ke[], ke_up)
            pe[] = push!(pe[], pe_up)
            energy[] = push!(energy[], e_up)
            reset_limits!(ax2)
            reset_limits!(ax3)
            reset_limits!(ax4)
            xlims!(ax2, 1, i+0.05)
            xlims!(ax3, 1, i+0.05)
            xlims!(ax4, 1, i+0.05)
            for i in eachindex(init_params.x_positions)
                drawx[i][] = [init_params.x_positions[i]] 
                drawy[i][] = [init_params.y_positions[i]] 
                drawz[i][] = [init_params.z_positions[i]] 
                new_point = Point3f0(init_params.x_positions[i], init_params.y_positions[i], init_params.z_positions[i])
                
                # Draws particle tails
                if length(line[i][]) > 300
                    line[i][] = push!(line[i][][end-300:end], new_point)
                else
                    line[i][] = push!(line[i][], new_point)
                end
                
            end
            println(init_params.Kinetic_Es[2],"   " ,init_params.Potential_Es[2],"   ", init_params.Potential_Es[2]+init_params.Kinetic_Es[2])
            ax1.title = "t = $(round(i*init_params.dt/86400, digits=4)) days"
            recordframe!(io)
        end
    end
end
plotme()



