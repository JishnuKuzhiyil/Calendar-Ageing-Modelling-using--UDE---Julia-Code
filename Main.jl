
using Pkg
Pkg.instantiate()



using DifferentialEquations,Plots,SparseArrays
using PreallocationTools,Statistics,MAT,LinearAlgebra
using ComponentArrays, Lux, DiffEqFlux, Plots.PlotMeasures,StableRNGs

#custom modules
include("model_parameters.jl") ; import .Para #contains electrochemical model parameters 
include("Cache_vectors.jl")    ; import .Cache_vectors  # contains Cache vectors to avoid heap allocation
include("Experiment.jl")       ; import .Experiment_func  # contains functions to define virtual experiments

SOC=5#Storage SOC
Temperature=45 #Storage Temperature
Model=("Physics","UDE")[2] #1 for Physics, 2 for UDE

dates,Exp_capacity,Exp_capacity_std = Experiment_func.get_capacity_data(Temperature,SOC)# Capacity data from RPT+RPTx
LAM_dates, LAM_mean, LAM_std = Experiment_func.get_LAM_data(Temperature,SOC) # LAM data from RPTx
Experiment = Experiment_func.Calendar_ageing_exp_from_dates(dates,SOC/100,Exp_capacity) #Experiment steps, implimented in Experiment_func Module

if Temperature==45 # Temperature dependent parameters

        k_SEI,α_SEI,D_SEI,U_SEI,β = 7.32e-16,0.5335,1.16e-21,0.4,4.526e-10 #Physics parameters
        κ₁ = 1.0;  κ₂=1.0 #UDE parameters

    elseif Temperature==25

        κ₁ = 0.2; κ₂=0.47 #UDE parameters
        k_SEI,α_SEI,D_SEI,U_SEI,β = 1.098e-16,0.5,1.856e-22,0.4,12.22e-10 # Physics parameters

    elseif Temperature==0

        κ₁ = 0.16 ; κ₂=0.3 # UDE parameters
        k_SEI,α_SEI,D_SEI,U_SEI,β = 2.928e-16,0.385,8.12e-23,0.4,13.85e-10 # Physics parameters

end 

#Defining Neural Networks using Lux.jl
sqr(x) = x.^2 #square activation fucntion
NN_eps =Lux.Chain(
    Lux.Dense(2, 10, tanh),
    Lux.Dense(10, 10, tanh),
    Lux.Dense(10, 1,sqr)
)
NN_SEI = Lux.Chain(
    Lux.Dense(2, 10, tanh),
    Lux.Dense(10, 10, tanh),
    Lux.Dense(10, 1,sqr)
)


rng = StableRNG(666) #random number generator
NN_p_SEI_, NN_st_SEI = Lux.setup(rng, NN_eps) #initialising neural networks with random parameters
NN_p_eps_, NN_st_eps = Lux.setup(rng, NN_eps) #initialising neural networks with random parameters
NN_p_SEI = ComponentArray(NN_p_SEI_) #converting to component array
NN_p_eps = ComponentArray(NN_p_eps_) #converting to component array
NN_p_SEI[1:end] = Para.NN_SEI_parameters # Changing ransom parameters to trained parameters
NN_p_eps[1:end] = Para.NN_eps_parameters # Changing ransom parameters to trained parameters


"""_______________________________________________Model Setup___________________________________________________________"""

    #Initial conditions, first 56 states are concentration, followed by SEIlength, negative porosity and anode active material fraction
    u0 = [Para.c_n_ini*ones(Para.N_FVM_particle); Para.c_p_ini*ones(Para.N_FVM_particle);Para.c_e_ref*ones(Para.N_FVM_electrode)]
    push!(u0,Para.L_SEI_ini) #SEI thickness
    push!(u0,Para.ϵ_neg_BOL) #Negative electrode porosity
    push!(u0,Para.ϵ_act_neg_BOL) #Negative electrode active material vol fraction
    push!(u0,0.0) #current at end (current become state function during voltage hold steps)

    dates_diff = diff(dates[1:end])
    Q_RPT=zeros(length(dates)) # Capacity will be stored in this array as RPTs are run
    global cycle_num=0 # RPT number
    T_exp=Temperature+273.15 # Storage Temperature in Kelvin

    #Parameters
    Exp_step_no=[1]
    Exp_step_transition=[0.0,0.0,0.0,0.0] #t_last,t_2ndlast,Voltage_last, Current_last

    #parameter set
    p=(
    
        #experiment related parameters
        Exp_step_no,Experiment,Exp_step_transition,

        #cache spaces to avoid heap allocation of arrays
        Cache_vectors.ϵ_vec_node,Cache_vectors.Brug_vec_node,
        Cache_vectors.grad_c_e_vector,Cache_vectors.brug_interface_vector,
        Cache_vectors.D_e_node,Cache_vectors.D_e_interface,
        Cache_vectors.arg_grad,Cache_vectors.grad_ce_flux,
        

        #SEI_parameters
        [α_SEI,k_SEI,D_SEI,U_SEI,β],

    )


    
"""___________________________________________Model Equations__________________________________________________________"""




function DAE_Equations!(du,u,p,t)


        @views c_n,c_p,c_e=u[1:10],u[11:20],u[21:56]        #Assign names to variables. use view macro to avoid array allocation
        @views dc_n,dc_p,dc_e=du[1:10],du[11:20],du[21:56]  #Assign names to state derivatives . use view macro to avoid array allocation
        L_SEI= u[57]                                        # SEI length state variable
        ϵ_neg = u[58]                                       # Negative electrode porosity state variable
        ϵ_act_neg= u[59]                                    # Negative electrode active material vol fraction state variable

        step=p[1][1]  #current experiment step number
        tlast=p[3][1]  # end time of last experimental step
        I = p[2][step][2](t-tlast) .- @view(u[end])  # Current at time t (the state u[end] is activated suring voltage hold)

        a_neg = Para.a_neg_BOL * (ϵ_act_neg/Para.ϵ_act_neg_BOL)  # Negative particle surface area density [m2/m3]
        F,R,T_ref = Para.F,Para.R,Para.T_ref


        """unpacking parameters___________________________________________________________"""

        ϵ_vec_node,Brug_vec_node = get_tmp(p[4],u),get_tmp(p[5],u)
        grad_c_e_vector,brug_interface_vector = get_tmp(p[6],u),get_tmp(p[7],u)
        D_e_node,D_e_interface = get_tmp(p[8],u),get_tmp(p[9],u)
        arg_grad,grad_ce_flux = get_tmp(p[10],u),get_tmp(p[11],u)

        α_SEI,k_SEI,D_SEI,U_SEI,β = p[12]

        """SEI calculations___________________________________________________________"""

        c_n_surf=1.5*c_n[end]-0.5*c_n[end-1]   #Linear extrapolation between last two nodes
        c_p_surf = 1.5*c_p[end]-0.5*c_p[end-1]  #Linear extrapolation between last two nodes
        OCP_n=Para.Negative_OCP(c_n_surf)     # Open circuit potential of anode
        OCP_p=Para.Positive_OCP(c_p_surf)     # Open circuit potential of cathode

        η_SEI_ohmic = I * Para.R_SEI *L_SEI     # Ohmic overpotential of SEI
        c_e_avg_n=mean(√, c_e[1:15])^2    # squared mean of root of electrolyte concentration
        j0_n=Para.Negative_exchange_current_density(c_e_avg_n,c_n_surf, Para.c_neg_max, T_ref) # Exchange current density of anode
        j_n = I/Para.A_electrode/Para.L_neg  # Current density of anode [A/m^3]
        η_SEI_reaction = 2*R*T_ref/F * asinh(j_n/j0_n/a_neg) # Reaction overpotential of SEI
        
        η_n_bar = OCP_n + η_SEI_reaction + η_SEI_ohmic # # Average overpotential of anode particle-SEI interface
        η_p_bar = OCP_p # Average overpotential of cathode particle-SEI interface  
        
      if Model=="Physics"

            SEI_exp_term = exp(-α_SEI*(F/R/T_exp)*(η_n_bar - U_SEI))  # Exponential term in SEI flux equation
            j_SEI =  -a_neg*F*Para.c0_SEI*k_SEI*SEI_exp_term/(1 + k_SEI*L_SEI*SEI_exp_term/D_SEI)  # SEI flux in [A/m³s]
            dL_SEI = -Para.Molar_weight_SEI/Para.rho_SEI/F/2/ Para.a_neg_BOL * j_SEI  # Rate of change of SEI thickness
            dϵ_neg = -1/a_neg*dL_SEI
            dϵ_act_neg = β*j_SEI

      elseif Model=="UDE"

            Ψ1(x) = leakyrelu(0.16-x,0.05)
            Ψ2(x) = exp(20(x-4.05))
            Ψ3(x) = exp(10.0 * leakyrelu(0.16-x,0.035)) 
            Ψ4(x) = x
            γ = 24*3600.0

            dL_SEI = 1e-9 * (1/γ)*NN_SEI([Ψ1(η_n_bar), Ψ2(η_p_bar)],NN_p_SEI,NN_st_SEI)[1][1]*κ₁/(1e9*L_SEI)
            dϵ_act_neg = -(1/100.0)*(1/γ )*NN_eps([Ψ3(η_n_bar), Ψ4(η_p_bar)],NN_p_eps,NN_st_eps)[1][1]*κ₂/(1+100.0*(Para.ϵ_act_neg_BOL-ϵ_act_neg))

            j_SEI =-a_neg * 2*Para.F*Para.rho_SEI/Para.Molar_weight_SEI* dL_SEI
            dϵ_neg = -1/a_neg*dL_SEI # Rate of change of negative electrode porosity

      end 
        du[57],du[58],du[59]=dL_SEI,dϵ_neg,dϵ_act_neg

        
        """Particle concentration equations__________________________________________"""

        mul!(dc_n,Para.NegParticle_StateMatrix.*0.6,c_n) #Particle diffusion: Anode, 0.6 is C rate correction factor for diffusivity.
        mul!(dc_p,Para.PosParticle_StateMatrix,c_p) #Particle diffusion: Cathode

        #Particle BCs
        j_p = -I/Para.A_electrode/Para.L_pos  # Current density of cathode [A/m^3]
        j_n = I/Para.A_electrode/Para.L_neg   # Current density of anode [A/m^3]

        N_p_surf = j_p/F/Para.a_pos_BOL # Flux at cathode particle surface [mol/m²s]
        N_n_surf = (j_n-j_SEI)/F/a_neg # Flux at anode particle surface [mol/m²s]

        dc_p[end]=dc_p[end].-Para.AbyV_outer_cntrl_vol_pos.*N_p_surf        #Cathode BC
        dc_n[end]=dc_n[end].-Para.AbyV_outer_cntrl_vol_neg.*N_n_surf        #Anode BC

    """Electrolyte concentration equations__________________________________________"""

        mul!(grad_c_e_vector,Para.Electrolyte_conc_interfaceGrad,c_e) # Gradient of electrolyte concentration
        ϵ_vec_node .= [ϵ_neg*ones(Para.N_FVM_neg_elec);Para.ϵ_sep_BOL*ones(Para.N_FVM_sep);Para.ϵ_pos_BOL*ones(Para.N_FVM_pos_elec)] # porosity at FVM nodes
        Brug_vec_node .= ϵ_vec_node.^1.5 # Brugman coefficient at FVM nodes
        mul!(brug_interface_vector,Para.Mean_of_node_at_interface,Brug_vec_node) # Bruggman at interface of control volumes
        D_e_node .= 8.794e-11 .* (c_e ./ 1000) .^ 2 .- 3.972e-10 .* (c_e ./ 1000) .+ 4.862e-10 # Diffusivity at FVM nodes
        mul!(D_e_interface,Para.Mean_of_node_at_interface,D_e_node) #Diffusivity at interface
        arg_grad .= D_e_interface .* grad_c_e_vector .* brug_interface_vector  # Flux at interfaces
        mul!(grad_ce_flux,Para.Electrolyte_divergence_matrix,arg_grad) # Divergence of flux at nodes
        dc_e .= (grad_ce_flux .+  Para.electrolyte_source*I)./ϵ_vec_node # Electrolyte source terms


        """Current and voltage algebraic expressions_________________________________________"""

        #Current dynamics
        if p[2][step][1]# during voltage hold
            V_terminal = Voltage_func(u,t,p)
            du[end]=V_terminal-p[2][step][6] #Voltage is constant
        else 
            du[end]=u[end]
        end


end 



function electrolyte_concentration_potential(ce;Implimentation = "Scot Moura's SPMeT")

      
    if Implimentation == "Ferran's SPMe-SR" 

        dx = [fill(85.2/15, 14); (85.2/15 + 2.0)/2; fill(2.0, 5); (2.0 + 75.6/15)/2; fill(75.6/15, 14)]

        dce_dx = [(ce[2] - ce[1]) / (dx[1] * 2);
                [(ce[i+1] - ce[i-1]) / (dx[i-1] + dx[i]) for i in 2:35];
                (ce[end] - ce[end-1]) /( dx[end] * 2)]
        
        norm_grad = dce_dx ./ ce
        integral = [0; cumsum(0.5 .* (norm_grad[1:end-1] .+ norm_grad[2:end]) .* dx)]
        integral[1] = 0.5 * norm_grad[1] * dx[1] / 2

    return coeff.*integral[1:15], coeff.*integral[22:36]

    elseif Implimentation == "Scot Moura's SPMeT"

        coeff = 2*Para.R*Para.T_ref/Para.F * (1-Para.t⁺) 
        return coeff.*log.(ce[1:15]./1000.0), coeff.*log.(ce[22:36]./1000.0)

    end 
    
    
end


function Voltage_func(u,t,p)


        @views c_n,c_p,c_e=u[1:10],u[11:20],u[21:56]        #Assign names to variables. use view macro to avoid array allocation
        L_SEI= u[57]                         # SEI length state variable
        ϵ_neg= u[58]                         # Negative electrode porosity state variable
        ϵ_act_neg= u[59]                  # Negative electrode active material vol fraction state variable

        step=p[1][1]  #current experiment step number
        tlast=p[3][1]  # end time of last experimental step
        I = p[2][step][2](t-tlast) .- @view(u[end])  # Current at time t
        

        a_neg = Para.a_neg_BOL * (ϵ_act_neg/Para.ϵ_act_neg_BOL)  # Negative particle surface area density [m2/m3]
        Brug_neg = ϵ_neg^1.5  # brugman coefficient for negative electrode

        OCP_pos = Para.Positive_OCP(1.5*c_p[end]-0.5*c_p[end-1])  # Open circuit potential of cathode
        OCP_neg = Para.Negative_OCP(1.5*c_n[end]-0.5*c_n[end-1])  # Open circuit potential of anode
        OCP = OCP_pos - OCP_neg  # Open circuit potential of cell


        Φe_conc_n = mean(electrolyte_concentration_potential(c_e;Implimentation = "Scot Moura's SPMeT" )[1])
        Φe_conc_p = mean(electrolyte_concentration_potential(c_e;Implimentation = "Scot Moura's SPMeT"  )[2])
        Φe_conc = Φe_conc_p - Φe_conc_n  # Concentration overpotential of electrolyte

        σe_neg_vec = Para.electrolyte_conductivity.(c_e[1:15])
        σe_pos_vec = Para.electrolyte_conductivity.(c_e[22:36])

        Φe_ohmic_e_n = mean((I*Para.L_neg/Para.A_electrode/Brug_neg/6)./σe_neg_vec)
        Φe_ohmic_e_p = -mean((I*Para.L_pos/Para.A_electrode/Para.Brug_pos/6)./σe_pos_vec)
        Φe_ohmic_e = Φe_ohmic_e_p - Φe_ohmic_e_n  # Ohmic overpotential of electrolyte

        Φs = -(I/Para.A_electrode )*Para.L_neg/Para.σ_neg - (I/Para.A_electrode )*Para.L_pos/Para.σ_pos  # Solid phase ohmic potential
        Φ_SEI = -I*L_SEI*Para.R_SEI  # SEI ohmic potential

        jn = I/Para.A_electrode/Para.L_neg  # Current density of anode [A/m^3]
        jp = -I/Para.A_electrode/Para.L_pos  # Current density of cathode [A/m^3]

        c_e_avg_n = mean(√, c_e[1:15])^2  # squared mean of root of electrolyte concentration
        j0_n = Para.Negative_exchange_current_density(c_e_avg_n,1.5*c_n[end]-0.5*c_n[end-1], Para.c_neg_max, Para.T_ref)  # Exchange current density of anode
        η_n_reaction = 2*Para.R*Para.T_ref/Para.F * asinh(jn/j0_n/a_neg)  # Overpotential of anode

        c_e_avg_p = mean(√, c_e[22:36])^2  # squared mean of root of electrolyte concentration
        j0_p = Para.Positive_exchange_current_density(c_e_avg_p,1.5*c_p[end]-0.5*c_p[end-1], Para.c_pos_max, Para.T_ref)  # Exchange current density of cathode
        η_p_reaction  = 2*Para.R*Para.T_ref/Para.F * asinh(jp/j0_p/Para.a_pos_BOL)  # Overpotential of cathode

        η_reaction = η_p_reaction  - η_n_reaction   # Reaction overpotential


        V = OCP + Φe_conc + Φe_ohmic_e + Φs + Φ_SEI + η_reaction  # Terminal voltage of cell

        return V
    
    
    end 


function condition(u, t, integrator)

    #experimental step number  is a parameter
    step=integrator.p[1][1]
    t_last=integrator.p[3][1]
    I = integrator.p[2][step][2](t-t_last) .- @view(u[end])

    experiment=integrator.p[2][:]
    exp_step=experiment[step]

    
    if exp_step[1]

        if exp_step[3]=="Time"

            return (t-t_last)-exp_step[4]

        elseif exp_step[3]=="Abs_Current"
            
            return abs(integrator.p[2][step][2](t-t_last) .- u[end]) .- exp_step[4]

        end

    else

        if exp_step[3]=="Time"

            return (t-t_last)-exp_step[4]

        elseif exp_step[3]=="Voltage"

            return Voltage_func(u,t,integrator.p)-exp_step[4]

        end

    end             



end 

function affect!(integrator)

    

    t_secondlast=integrator.p[3][1]
    t_last=integrator.t  
    step_nu=integrator.p[1][1]
    I_last=integrator.p[2][step_nu][2](t_last - t_secondlast)-integrator.u[end]

    if integrator.p[2][step_nu][6] == "RPT-Capacity"
        global cycle_num=cycle_num+1
        println("Running RPT number : $(cycle_num) out of $(length(dates)) RPTs. ")
        global Q_RPT[cycle_num] = (t_last-t_secondlast).*I_last./3600
    end

    if step_nu != length(integrator.p[2])
    integrator.opts.dtmax=integrator.p[2][step_nu+1][7] #update maximum time step for next experimental step
    end 
    
    if step_nu == length(integrator.p[2])

            terminate!(integrator)

    else

        if integrator.p[2][step_nu][1]

            integrator.u[end]=0.0
            nxt_stp_duration=integrator.p[2][step_nu+1][5]
            integrator.p[1][1]=step_nu.+1
            integrator.p[3].=[t_last,t_secondlast,0.0,I_last]
            reinit!(integrator,integrator.u,t0=t_last,tf=t_last+nxt_stp_duration,erase_sol=false)

        elseif integrator.p[2][step_nu+1][1]
        
            
            nxt_stp_duration=integrator.p[2][step_nu+1][5]
            integrator.p[3].=[t_last,t_secondlast,0.0,I_last]
            integrator.p[1][1]=step_nu.+1
            #integrator.u[end]=I_last
        
            reinit!(integrator,integrator.u,t0=t_last,tf=t_last+nxt_stp_duration,erase_sol=false)

        else
        
            integrator.p[1][1]=step_nu.+1
            integrator.p[3].=[t_last,t_secondlast,0.0,I_last]
            nxt_stp_duration=integrator.p[2][step_nu+1][5]

            reinit!(integrator,integrator.u,t0=t_last,tf=t_last+nxt_stp_duration,erase_sol=false)
        end
    end
    
end

cb=ContinuousCallback(condition,affect!)


"""______________________________________________solving the model___________________________________________________"""




Inertia_matrix=spdiagm([ones(length(u0)-1);0]) # singular mass matrix
tspan=(0.0,3.5*3600.0) #time span larger than time to first event
DAE_system=ODEFunction(DAE_Equations!,mass_matrix=Inertia_matrix)

prob=ODEProblem(DAE_system,u0,tspan,p)
sol=solve(prob,Rosenbrock23(),callback=cb,dtmax=300.0) 

"""______________________________________________Post processing___________________________________________________"""

sol_ϵ_act_neg(day) = sol(day*24*3600)[59]
LAM_sim = 100.0 .- 100.0 .*[sol_ϵ_act_neg(day) for day in LAM_dates]./Para.ϵ_act_neg_BOL
Q_RPT_norm = Q_RPT./Q_RPT[1] .*100.0
Exp_capacity_norm = Exp_capacity./Exp_capacity[1] .*100.0
Exp_capacity_std_norm = Exp_capacity_std./Exp_capacity[1] .*100.0


#plotting the capacity results
f=16
P =plot()
scatter!(P,dates,Q_RPT_norm, label=Model,markersize=6,markercolor=:green,marker=:diamond)
scatter!(P,dates,Exp_capacity_norm, yerr=Exp_capacity_std_norm, label="Experiment",markersize=6,marker=:circle,markercolor=:black,legend=:bottomleft)
plot!(P,framestyle=:box,xlabel="Days of storage",ylabel="Relative Capacity (%)")
plot!(P,xlabelfontsize=f,ylabelfontsize=f,xtickfontsize=f-4,ytickfontsize=f-4,legendfontsize=f-4,titlefontsize=f)  

#plotting the LAM results
P2=plot()
scatter!(P2,LAM_dates,LAM_sim, label=Model,markersize=6,markercolor=:green,marker=:diamond)
scatter!(P2,LAM_dates,LAM_mean, yerr=LAM_std',label="Experiment",markersize=6,marker=:circle,markercolor=:black,legend=:topleft)
plot!(P2,framestyle=:box,xlabel="Days of storage",ylabel="LAM (%)")
plot!(P2,xlabelfontsize=f,ylabelfontsize=f,xtickfontsize=f-4,ytickfontsize=f-4,legendfontsize=f-4,titlefontsize=f)

RMSE_cap = sqrt(mean((Exp_capacity_norm-Q_RPT_norm).^2))
RMSE_LAM = sqrt(mean((LAM_mean-LAM_sim).^2))
println("RMSE for capacity is $(RMSE_cap) and for LAM is $(RMSE_LAM)")


plot(P,P2,layout=(2,1),size=(500,600))