using ITensors
using PyCall
using NPZ

let
    num=151
    N_point=Float64[]
    E_point=Float64[]

    N = 201
    sites = siteinds("S=1/2",N)

    hz=0.1
    m=0.1
    a=1
    os_N = OpSum()
    for j=1:N-1
        os_N += -2,"Sz",j,"Sz",j+1
    end
    NN = MPO(os_N,sites)

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

    OE = MPO(os_E,sites)

    for k=1:num
        mu=0.118+0.002*k
        os = OpSum()

        for j=1:N
            os += 2*hz,"Sz",j
        end

        for j=1:N-1
            os += 2*(mu-m*((-1)^j)),"Sz",j,"Sz",j+1
        #os += 4*mu,"Sz",j,"Sz",j+1
        end

        for j=1:N-2
            os += -2/a,"Sz",j,"Sx",j+1,"Sz",j+2
            os += 1/2/a,"Sx",j+1
        end

        H = MPO(os,sites)

        nsweeps = 300 # number of sweeps is 5
        maxdim = [10,20,100,100,200,200,500,1000,1500] # gradually increase states kept
        cutoff = [1E-10]

        psi = randomMPS(sites,1024)

        energy,psi = dmrg(H,psi; nsweeps, maxdim, cutoff)

        N_val=real(inner(psi',NN,psi))
        push!(N_point,N_val)

        E_val=real(inner(psi',OE,psi))
        push!(E_point,E_val)
    end

    file_name="Data_E_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)
    npzwrite(file_name, E_point)
    file_name="Data_N_h_"*string(hz)*"_a_"*string(a)*"_L"*string(N)
    npzwrite(file_name, N_point)
end

#Aver_h_corr=h_data_point 



