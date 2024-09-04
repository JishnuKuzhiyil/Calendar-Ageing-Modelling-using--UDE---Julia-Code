# model_parameters.jl

module Para 

using LinearAlgebra,SparseArrays,PreallocationTools,Statistics
# Global parameters
const T_ref = 298.15  # Refernce Temperature[K]
const F = 96485.3329  # Faraday constant [C/mol]
const R = 8.314       # Ideal gas constant [J/(mol K)]



# Electrode parameters
const L_pos = 75.6e-6  # Positive electrode thickness [m]
const L_neg = 85.2e-6  # Negative electrode thickness [m]
const L_sep = 12.0e-6  # Separator thickness [m]
const R_pos = 5.22e-6  # Positive particle radius [m]
const R_neg = 5.86e-6  # Negative particle radius [m]
const A_electrode = 0.1 # Area of electrode [m^2]

const ϵ_pos_BOL = 0.335  # Positive electrode porosity
const ϵ_neg_BOL = 0.25  # Negative electrode porosity
const ϵ_sep_BOL = 0.47   # Separator porosity
const Brug_pos = ϵ_pos_BOL ^1.5  # Bruggeman coefficient in the positive electrode
const Brug_sep = ϵ_sep_BOL ^1.5  # Bruggeman coefficient

const a_neg_BOL = 3.84e5 # Negative particle surface area density [m^2/m^3]
const a_pos_BOL = 3.82e5 # Positive particle surface area density [m^2/m^3]
const ϵ_act_neg_BOL = 0.75  # Active material volume fraction in the negative electrode
const ϵ_act_pos_BOL = 0.665  # Active material volume fraction in the positive electrode

const N_FVM_pos_elec = 15  # Number of FVM elements in the positive electrode
const N_FVM_neg_elec = 15  # Number of FVM elements  in the negative electrode
const N_FVM_sep = 6  # Number of FVM elements in the separator
const N_FVM_electrode = N_FVM_pos_elec + N_FVM_neg_elec+ N_FVM_sep # Total number of FVM elements in the electrode
const N_FVM_particle =10 # Number of FVM elements in the particle (anod eand cathode)
const AbyV_outer_cntrl_vol_neg = 1.889097389e6 # Area to volume ratio of the last control volume in the negative particle
const AbyV_outer_cntrl_vol_pos = 2.120710862e6 # Area to volume ratio of the last control volume in the positive particle

const c_neg_max  =33133.0   # Maximum concentration of lithium in the negative particle [mol/m^3]
const c_pos_max  =  63104.0 # Maximum concentration of lithium in the positive particle [mol/m^3]
const c_e_ref    = 1000.0   # Reference concentration of lithium in the electrolyte [mol/m^3]

const t⁺ = 0.2594 # Transference number of lithium ions in the electrolyte
σ_neg = 0.18 # Conductivity of the negative electrode [S/m]
σ_pos = 215 # Conductivity of the positive electrode [S/m]

SOC_ini = 1.0         # Initial state of charge
θn_SOC0,θn_SOC100 = 0.02637,0.910612
θp_SOC0,θp_SOC100 = 0.8539736,0.263848
c_n_ini = c_neg_max *(θn_SOC0+(θn_SOC100-θn_SOC0)*SOC_ini)    #mol/m^3
c_p_ini = c_pos_max*(θp_SOC0+(θp_SOC100-θp_SOC0)*SOC_ini)   #mol/m^3
L_SEI_ini = 1e-9 # Initial thickness of SEI layer [m]

const R_SEI = 6.1e4 # SEI resistance [Ohm m^2]
const c0_SEI = 4541.0 # solvent concentration [mol/m^3]
const Molar_weight_SEI = 0.162 # Molar weight of solvent in SEI [kg/mol]
const rho_SEI = 1690 # Density of SEI [kg/m^3]



#Positive particle state matrix (Tridiagonal sparse matrix)(D_p=4e-15 m^2/s , 10 FVM elements, R_p=5.22e-6 m )
PosParticle_LDiag = [0.006291326148638871, 0.009271428008520441, 0.010712258036871591, 0.011551287354877928, 0.012098704131997828, 0.012483576294936972, 0.012768786207000784, 0.012988544306867347, 0.013163032938296088]
PosParticle_Diag=[-0.044039283040472096, -0.03145663074319435, -0.030132141027691434, -0.02975627232464331, -0.02960017384687469, -0.029520838082074702, -0.029475110696378964, -0.029446384518185483, -0.029427170695246332, -0.013163032938296088]   
PosParticle_UDiag= [0.044039283040472096, 0.025165304594555484, 0.020860713019170994, 0.019044014287771718, 0.018048886491996763, 0.01742213395007687, 0.01699153440144199, 0.016677598311184698, 0.016438626388378987]
const PosParticle_StateMatrix=Tridiagonal(PosParticle_LDiag,PosParticle_Diag,PosParticle_UDiag)

#Negative particle state matrix (Tridiagonal sparse matrix)(D_n=3.3e-14 m^2/s , 10 FVM elements, R_p=5.86e-6 m )
const NegParticle_StateMatrix = ( 6.546357558037951) .* PosParticle_StateMatrix


dx_neg = L_neg/N_FVM_neg_elec
dx_pos = L_pos/N_FVM_pos_elec
dx_sep = L_sep/N_FVM_sep
dx_neg_sep = (dx_neg + dx_sep)/2
dx_pos_sep = (dx_pos + dx_sep)/2

e_vec = 1 ./[dx_neg.*ones(N_FVM_neg_elec-1);dx_neg_sep;dx_sep.*ones(N_FVM_sep-1);dx_pos_sep;dx_pos.*ones(N_FVM_pos_elec-1)]
const Electrolyte_conc_interfaceGrad=spdiagm(N_FVM_electrode-1,N_FVM_electrode,0=>-e_vec,1=>e_vec)

left_weight = 0.5*ones(N_FVM_electrode-1)
right_weight = 0.5*ones(N_FVM_electrode-1)
left_weight[N_FVM_neg_elec],right_weight[N_FVM_neg_elec] = dx_neg/(dx_neg+dx_sep),dx_sep/(dx_neg+dx_sep)
left_weight[N_FVM_neg_elec+N_FVM_sep],right_weight[N_FVM_neg_elec+N_FVM_sep] = dx_sep/(dx_sep+dx_pos),dx_pos/(dx_sep+dx_pos)
const Mean_of_node_at_interface = spdiagm(N_FVM_electrode-1,N_FVM_electrode,0=>left_weight,1=>right_weight)

left_weight1 = 1 ./[dx_neg.*ones(N_FVM_neg_elec-1);dx_sep*ones(N_FVM_sep);dx_pos*ones(N_FVM_pos_elec)]
right_weight1 = 1 ./[dx_neg*ones(N_FVM_neg_elec);dx_sep*ones(N_FVM_sep);dx_pos*ones(N_FVM_pos_elec-1)]
const Electrolyte_divergence_matrix = spdiagm(N_FVM_electrode,N_FVM_electrode-1,0=>right_weight1,-1=>-left_weight1)


Source_neg = (1-t⁺)/(F*L_neg*A_electrode)
Source_sep = 0.0
Source_pos = -(1-t⁺)/(F*A_electrode*L_pos)
const electrolyte_source = [Source_neg*ones(N_FVM_neg_elec);Source_sep*ones(N_FVM_sep);Source_pos*ones(N_FVM_pos_elec)]




function Positive_OCP(surf_conc)
        
    sto=surf_conc/c_pos_max

    OCP=(
        -0.8090 * sto
        + 4.4875
        - 0.0428 * tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * tanh(15.9308 * (sto - 0.3120))
    )
    return OCP
end 


function Negative_OCP(surf_conc)

    sto=surf_conc/c_neg_max

    OCP = (
        1.9793 * exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 *tanh(29.8538 * (sto - 0.1234))
        - 0.04478 *tanh(14.9159 * (sto - 0.2769))
        - 0.0205 *tanh(30.4444 * (sto - 0.6103))
        )
    return OCP
end 

function Negative_exchange_current_density(c_e, c_s_surf, c_s_max, T)

    #T=298.15
    m_ref = 2.1e-5 #6.48e-7 # # (A/m2)(m3/mol)^1.5 - includes ref concentrations
    E_r = 35000
    arrhenius = exp(E_r / 8.314 * (1 / 298.15 - 1 ./ T))
   
    if c_s_surf<0 || (c_s_max - c_s_surf) <0
        return  -(m_ref * arrhenius * c_e^0.5 * abs(c_s_surf)^0.5 * abs((c_s_max - c_s_surf)) ^ 0.5)
    end

    return  m_ref * arrhenius * c_e^0.5 * c_s_surf^0.5 * (c_s_max - c_s_surf) ^ 0.5
end

function Positive_exchange_current_density(c_e, c_s_surf, c_s_max, T)

    m_ref = 3.42e-6  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    E_r = 17800
    arrhenius = exp(E_r / 8.314 * (1 / 298.15 - 1 ./ T))

    if c_s_surf<0 || (c_s_max - c_s_surf) <0
        return  -(m_ref * arrhenius * c_e^0.5 * abs(c_s_surf)^0.5 * abs((c_s_max - c_s_surf)) ^ 0.5)
    end
    return m_ref * arrhenius * c_e^0.5 * c_s_surf^0.5 * (c_s_max - c_s_surf)^0.5
    
end 


function electrolyte_conductivity(ce)

    return 1.297e-10*ce^3 - 7.937e-5*ce^1.5 + 3.329*ce 
end 

#Nural network parameters: [1:20](layer1_weights), [21:30](layer1_bias),[31:130](layer2_weights),[131:140](layer2_bias),[141:150](layer3_weights),[151:151](layer3_bias)
NN_SEI_parameters= Float32[0.22198194 1.7313771 -0.47964984 -0.94933736 -2.7753005 0.43130958 1.8931949 -0.16180421 0.23406927 -1.280899 0.9603463 -0.3343974 -0.22627877 -0.021987284 -0.09193185 0.24473588 -0.37687436 0.1136438 -0.25508964 0.05237186 0.10539822 -0.54830015 0.46892494 -0.23345149 0.71488595 -0.3999595 -0.72287196 -0.85755694 0.30351835 0.43303996 0.15370524 0.08290501 0.006461135 -0.109199785 0.59386027 0.26953092 0.39379746 -0.6425541 0.23654316 -0.09581907 -0.5741323 0.078115724 -0.061593775 -0.96632975 0.98748064 0.17144574 -0.39169404 -0.23438638 0.09842143 -0.32691753 -0.36364508 0.7148075 0.50605434 0.27635336 -0.3959376 -0.4774318 0.5827349 0.27437353 -0.040344536 0.58928263 -0.1746948 -0.52379286 -0.44447783 0.5989246 -0.9633133 -0.4810085 -0.48930395 0.00487655 -0.13752241 -0.26354086 0.561412 0.28603768 -0.38607842 1.378465 -1.7087533 -0.1442259 0.16642553 0.29893982 -1.0923637 0.5724012 0.39859617 -0.34349135 0.09347034 -0.21653727 -0.06512126 0.16055351 -0.3084221 -0.537241 0.31552264 0.124487475 -0.14778017 -0.5692065 0.7380955 -0.7110882 0.98488855 0.71164906 -0.51833385 -0.21919677 0.7862283 -0.046964478 -0.5197486 -0.32457197 0.47228238 0.31475493 -0.19231087 0.08660056 -0.21617216 0.066772774 0.024403637 -0.48482388 -0.13798611 0.23538515 0.17841995 -0.045454536 -0.09284309 0.3249392 0.42155406 0.17649825 -0.18938798 0.36196098 -0.017211057 0.31437078 -0.5335321 0.5316772 -0.80758053 0.033115022 0.4355313 0.34980568 0.14076328 0.5997589 0.45379102 0.6321895 -0.31913444 0.20082353 -0.20017406 -0.23872702 0.47214296 0.21003063 -0.4602013 0.2675331 0.25444305 0.7703698 -0.22343019 0.8091879 -1.3121551 -0.006489412 0.4452825 0.2620436 -0.52950835 0.09093349 0.36747158]'
NN_eps_parameters = Float32[0.41589388 0.21131223 0.3501634 -0.5580254 -0.51451504 -0.23254189 0.23103195 -0.68229944 -0.12394388 -0.15580808 0.6555258 -0.39541858 -0.261224 0.335405 0.18642956 0.08558504 -0.44232872 0.13382743 -0.38493416 0.10446707 0.0036284912 0.021559887 -0.05785473 0.026371736 0.0069512906 -0.004343542 0.021920219 0.016172603 -0.01078714 0.002077939 0.21953993 -0.0015278115 -0.08248645 0.11642474 0.35048732 0.31206915 0.26818216 -0.48677036 0.113083735 -0.22133832 -0.26101458 0.34394556 -0.3717673 -0.4376884 0.39596704 -0.17354557 -0.1358698 -0.082761526 -0.17099553 -0.06088108 -0.4007425 0.57222575 0.5327186 0.2932601 -0.37896812 -0.41873065 0.41805935 0.23109409 0.09853279 0.45613667 -0.31400868 -0.45835614 -0.27062884 0.09692758 -0.36758465 -0.37766987 -0.41066507 -0.094406575 -0.059593957 -0.23601688 -0.08995362 -0.15050551 0.31199068 0.2783065 -0.49815854 0.50001574 -0.25091863 -0.15735427 -0.36507958 0.1217839 0.54287136 -0.07994394 -0.008011857 -0.21759592 -0.07607581 0.02198761 -0.058512967 -0.43714264 0.16607076 0.33015814 0.16359681 -0.31658426 0.4050467 -0.09695607 0.2845252 0.37399957 -0.25280178 -0.032470822 0.38820615 0.21326749 -0.5402334 -0.1387367 0.4924812 0.13744003 -2.0515858f-5 0.07068752 -0.025440728 0.07883726 -0.12069203 -0.34648338 -0.16533539 0.17979297 0.18337324 -0.01653489 -0.14540067 0.23966268 0.39725572 0.13230625 -0.10447101 0.40599167 -0.37351787 -0.04588231 -0.20412815 0.07582083 -0.31657133 0.36370456 0.100928664 0.12419752 0.4537495 0.2842876 0.0056735915 0.0037727044 0.0035641494 -0.0023810535 0.007740005 0.0042070537 -0.014876971 -0.00085576397 0.0034569406 -0.0035772915 0.43332016 0.6297577 -0.28968772 0.6990095 -0.62820303 -0.030963605 0.3366027 0.3100371 -0.21821636 0.051840153 0.031859126]'



end 