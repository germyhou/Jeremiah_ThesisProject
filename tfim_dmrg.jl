using ITensors, ITensorMPS, Plots
using DelimitedFiles 
using CSV
using DataFrames

let
    L=100
    h=0.5
    hb=sqrt(1-h)
    J=1
    sites = siteinds("S=1/2", L)

    function create_os(L, h, hb, J)
        os = OpSum()
        for n = 1:L-1
            os += -4*J,"Sz",n,"Sz",n+1
        end
        for n = 1:L
            os += -2*h,"Sx",n
        end
        os += -2*hb,"Sz",1
    end

    os = create_os(L, h, hb, J)
    H = MPO(os, sites)


    state = [isodd(n) ? "Up" : "Dn" for n in 1:L]
    psi0_init = productMPS(sites, state)

    nsweeps = 10
    maxdim = [10, 20, 50, 100, 100, 200, 200, 200, 400, 400]
    cutoff = 1E-8
    noise = [1E-6, 1E-6, 1E-7, 1E-8, 0, 0, 0, 0]

    energy, psi = dmrg(H, psi0_init; nsweeps=nsweeps, maxdim=maxdim, cutoff=cutoff, noise=noise)
    println("\nGround state energy = ", energy)
    println("\nTotal Sz of Ground State = ", sum(expect(psi, "Sz")))
    println("\nOverlap with the initial state = ", inner(psi0_init, psi))

    tau = 0.05
    cutoff = 1E-8
    ttotal= 100

    gates = ITensor[]
    for j in 1:(L-1)
      s1 = sites[j]
      s2 = sites[j + 1]
      hj =
        -4*J * op("Sz", s1) * op("Sz", s2)
      Gj = exp(-im * tau / 2 * hj)
      push!(gates, Gj)
    end
    for j in 1:L
        hj = -2*h * op("Sx",sites[j])
        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
    end

    h_last = -2*hb * op("Sz",sites[1])
    G_last = exp(-im * tau / 2 * h_last)
    push!(gates, G_last)
    append!(gates, reverse(gates))

    # <ψ|σz1|ψ>
    exp_σz1 = 2*expect(psi, "Sz", sites=1)

    psi0 = copy(psi)

    sop = AutoMPO()
    sop += 2, "Sz", 1
    Operator = MPO(sop, sites)
    finalize(sop)
    finalpsi = apply(Operator, psi)
    sp_psi = noprime(finalpsi)
    finalize(finalpsi)
    finalize(psi)
    println("Computing dynamics")

    corrdata = []
    data = []

    # Compute and print <Sz> at each time step
    # then apply the gates to go to the next time
    for t in 0.0:tau:ttotal
        # Sz = expect(sp_psi, "Sz", sites=1)
        # println("$t $Sz")

        corr_firstterm = exp(-im*energy*t)*inner(sp_psi,finalpsi)
        corr = corr_firstterm - exp_σz1^2
        println("$t $corr")
        push!(corrdata,corr)
        push!(data, (t, corr))

        t≈ttotal && break

        sp_psi = apply(gates, sp_psi; cutoff)
        normalize!(sp_psi)
    end

    # Save the time and correlation data to a text file
    # writedlm("time_correlation_data.txt", data, ' ')  # Write data with space delimiter

    time = 0.0:tau:ttotal

    y_axis = log.(abs.(corrdata))
    x_axis = log.(time)

    df = DataFrame(time = x_axis, data = y_axis)
    CSV.write("output.csv", df)
    
end