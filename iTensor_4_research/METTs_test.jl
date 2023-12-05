using ITensors
using PyCall
using NPZ

function avg_err(v::Vector)
    N = length(v)
    avg = v[1] / N
    avg2 = v[1]^2 / N
    for j in 2:N
      avg += v[j] / N
      avg2 += v[j]^2 / N
    end
    return avg, √((avg2 - avg^2) / N)
end

let
    #function main(; N=10, cutoff=1E-8, δτ=0.1, beta=2.0, NMETTS=3000, Nwarm=10)
    
    mu_step=parse(Int64, ARGS[1])
    N=101
    cutoff=1E-9
    δτ=0.05
    beta=30.0
    NMETTS=500
    Nwarm=10
    

    hz=0.1
    m=0.1
    mu=0.12+mu_step*0.01
    a=1
    # Make an array of 'site' indices
    s = siteinds("S=1/2", N)
    gates = ITensor[]

    for j in 1:(N)
        s1 = s[j]
        hi = op("Sz", s1)
        Gi = exp(-hz*δτ* hi)
        
        push!(gates, Gi)
    end

    for j in 1:(N - 2)
        s1 = s[j]
        s2 = s[j + 1]
        s3 = s[j + 2]
        hi = op("Id", s1)*op("Sx", s2)* op("Id", s3)+ -4*op("Sz", s1) * op("Sx", s2) * op("Sz", s3)

        Gi = exp(-δτ/a/4 * hi)
        
        push!(gates, Gi)
    end

    for j in 1:(N - 1)
        s1 = s[j]
        s2 = s[j + 1]

        hj = op("Sz", s1) * op("Sz", s2)
        Gj = exp(-1*(mu-m*(-1)^j)*δτ* hj)
        push!(gates, Gj)
    end


    # Include gates in reverse order too
    # (N,N-1),(N-1,N-2),...
    append!(gates, reverse(gates))

    # Initialize psi to be a product state (alternating up and down)
    #psi = MPS(s, n -> isodd(n) ? "Up" : "D
    # Make gates (1,2),(2,3),(3,4),...
    #gates = ops([("expτSS", (n, n + 1), (τ=-δτ / 2,)) for n in 1:(N - 1)], s)
    # Include gates in reverse order to complete Trotter formula
    #append!(gates, reverse(gates))

    # Make y-rotation gates to use in METTS collapses
    Rz_gates = ops([("Rz", n, (θ=-π / 2,)) for n in 1:N], s)
    Ry_gates = ops([("Ry", n, (θ=-π / 2,)) for n in 1:N], s)

    # Arbitrary initial state
    psi = randomMPS(s)

    # Make H for measuring the energy
    os = OpSum()

    for j=1:N
        os += 2*hz,"Sz",j
    end

    for j=1:N-1
        os += 2*(mu-m*((-1)^j)),"Sz",j,"Sz",j+1
    end

    for j=1:N-2
        os += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
        os += 1/2/a,"Sx",j+1
    end

    H = MPO(os,s)
    #H = MPO(terms, s)

    os_N = OpSum()
    for j=1:N-1
        os_N += -2,"Sz",j,"Sz",j+1
    end
    NN = MPO(os_N,s)

    NNN=apply(NN, NN; cutoff)

    os_E = OpSum()

    for j=1:N
        os_E += 2*hz,"Sz",j
    end

    for j=1:N-1
        os_E += 2*(-m*((-1)^j)),"Sz",j,"Sz",j+1
    #os += 4*mu,"Sz",j,"Sz",j+1
    end

    for j=1:N-2
        os_E += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
        os_E += 1/2/a,"Sx",j+1
    end

    OE = MPO(os_E,s)

    # Make τ_range and check δτ is commensurate
    τ_range = δτ:δτ:(beta / 2)
    if norm(length(τ_range) * δτ - beta / 2) > 1E-10
        error("Time step δτ=$δτ not commensurate with beta/2=$(beta/2)")
    end

    energies = Float64[]
    Num_lst = Float64[]
    N_s_lst = Float64[]

    for step in 1:(Nwarm + NMETTS)
        if step <= Nwarm
            println("Making warmup METTS number $step")
        else
            println("Making METTS number $(step-Nwarm)")
        end

    # Do the time evolution by applying the gates
        for τ in τ_range
            psi = apply(gates, psi; cutoff)
            normalize!(psi)
        end

    # Measure properties after >= Nwarm 
    # METTS have been made
        if step > Nwarm
            energy = inner(psi', OE, psi)
            Nums = inner(psi', NN, psi)
            Num_sq = inner(psi', NNN, psi)
            push!(energies, real(energy))
            push!(Num_lst, real(Nums))
            push!(N_s_lst, real(Num_sq))
            
            println("  Energy of METTS ", energy)
            a_E, err_E = avg_err(energies)
            println("Estimated Energy = ",a_E)
            println("  Number of METTS ", real(Nums))
            a_N, err_N = avg_err(Num_lst)
            println("Estimated number = ",a_N)
            println(" Squared Number of METTS ", real(Num_sq))
            a_NN, err_NN = avg_err(N_s_lst)
            println("Estimated square number = ",a_NN)
        end

    # Measure in X or Z basis on alternating steps
        #if step % 2 == 1
        psi = apply(Rz_gates, psi)
        psi = apply(Ry_gates, psi)
        #  samp = sample!(psi)
        #  new_state = [samp[j] == 1 ? "X+" : "X-" for j in 1:N]
        #else
        #  samp = sample!(psi)
        #  new_state = [samp[j] == 1 ? "Z+" : "Z-" for j in 1:N]
        #end
        samp = sample!(psi)
        new_state = [samp[j] == 1 ? "Y+" : "Y-" for j in 1:N]
        psi = productMPS(ComplexF64, s, new_state)
    end
    mu=round(mu, digits=3)
    file_name="Data_E_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)*"_mu_"*string(mu)*"_500"
    npzwrite(file_name, energies)
    file_name="Data_N_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)*"_mu_"*string(mu)*"_500"
    npzwrite(file_name, Num_lst)
    file_name="Data_Ns_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)*"_mu_"*string(mu)*"_500"
    npzwrite(file_name, N_s_lst)

end

#Aver_h_corr=h_data_point 



