using GLMakie

const G = 6.67e-11
Base.@kwdef mutable struct Params
    x_positions::Array{Float64}
    y_positions::Array{Float64}
    masses::Array{Float64}
    x_velocities::Array{Float64}
    y_velocities::Array{Float64}
    x_accelerations::Array{Float64}
    y_accelerations::Array{Float64}
end
nbody = 80
min_lim = -5e3
max_lim = 5e3
axis_lim_min = -1e4
axis_lim_max = 1e4
pos_x = rand(min_lim:max_lim, nbody)
pos_y = rand(min_lim:max_lim, nbody)
masses = rand(1:0.1:2, nbody)
x_velocities = zeros(nbody)
y_velocities = zeros(nbody)
x_accelerations = zeros(nbody)
y_accelerations = zeros(nbody)

# y_velocities = rand(-1:1, nbody)
# x_velocities = rand(-1:1, nbody)
init_params = Params(pos_x, pos_y, masses, x_velocities, y_velocities, x_accelerations, y_accelerations)

function get_acc(par::Params)
    for i in eachindex(par.x_positions)
        A_x, A_y = 0,0
        for p in eachindex(par.x_positions)
            if p != i
                #println(par.x_positions[i] - par.x_positions[p])
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 +50^2)
                #println("Radius: ", r)
                A_x += -masses[p]*(par.x_positions[i]-par.x_positions[p])/r^3
                A_y += -masses[p]*(par.y_positions[i]-par.y_positions[p])/r^3
                par.x_accelerations[i] = A_x
                par.y_accelerations[i] = A_y
            end
        end
    end
end

function init_velocities(par::Params)
    for i in eachindex(par.x_velocities)
        dt = 0.1
        for p in eachindex(par.x_velocities)
            if p!=i
                r = sqrt.( (par.x_positions[i] - par.x_positions[p])^2 + (par.y_positions[i] - par.y_positions[p])^2 )
                Vx = par.x_velocities[i] - dt*0.5G*masses[p]*(par.x_positions[i]-par.x_positions[p])/r^1.5
                Vy = par.y_velocities[i] - dt*0.5G*masses[p]*(par.y_positions[i]-par.y_positions[p])/r^1.5
                par.x_velocities[i] = Vx
                par.y_velocities[i] = Vy
                
            end
        end
    end
end


function leapfrog(par::Params)
    dt=20
    for i in eachindex(par.x_positions)
        Vx_update = par.x_velocities[i] + dt*par.x_accelerations[i]
        Vy_update = par.y_velocities[i] + dt*par.y_accelerations[i]
        X_pos_update = par.x_positions[i] + dt*Vx_update
        Y_pos_update = par.y_positions[i] + dt*Vy_update
        par.x_velocities[i] = Vx_update
        par.y_velocities[i] = Vy_update
        par.x_positions[i] = X_pos_update
        par.y_positions[i] = Y_pos_update
    end
end

fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (1200, 800))
ax1 = Axis(fig[1, 1])
limits!(ax1, axis_lim_min, axis_lim_max, axis_lim_min, axis_lim_max)


drawx,drawy = [], []
line = []
for i in eachindex(init_params.x_positions)
    append!(drawx, [Node([init_params.x_positions[i]])])
    append!(drawy, [Node([init_params.y_positions[i]])])
    append!(line, [Node(Point2f0[(init_params.x_positions[i][], init_params.y_positions[i][])])])
    lines!(line[i])
    scatter!(drawx[i], drawy[i], markersize = init_params.masses[i]*5)
end


record(fig, "test.mp4") do io
    for i in 1:100000
        get_acc(init_params)
        leapfrog(init_params)
        X_max = maximum(init_params.x_positions)
        X_min = minimum(init_params.x_positions)
        Y_max = maximum(init_params.y_positions)
        Y_min = minimum(init_params.y_positions)

        for i in eachindex(init_params.x_positions)
            drawx[i][] = [init_params.x_positions[i]] 
            drawy[i][] = [init_params.y_positions[i]] 
            new_point = Point2f0(init_params.x_positions[i], init_params.y_positions[i])
            if length(line[i][]) > 20
                line[i][] = push!(line[i][][end-20:end], new_point)
            else
                line[i][] = push!(line[i][], new_point)
            end
        end
        recordframe!(io)
    end
end
