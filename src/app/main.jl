include("base.jl")

function main()
    sys = RosslerSystem(0.1, 0.1, 12.0, 0.00005)
    backend = Tsit5()
    sol = integrate!(sys; tspan=(0.0, 20000.0), saveat=0.01)
end
