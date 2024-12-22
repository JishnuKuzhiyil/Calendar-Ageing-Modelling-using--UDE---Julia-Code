module Cache_vectors

    using LinearAlgebra,SparseArrays,PreallocationTools,Statistics,Revise
    includet("Model_parameters.jl")
    import ..Para

    N_elec = Para.N_FVM_electrode

    ϵ_vec_node = DiffCache(zeros(Float64,N_elec)) # Electrode porosity
    Brug_vec_node = DiffCache(zeros(Float64,N_elec)) # brugeman porosity

    grad_c_e_vector = DiffCache(zeros(Float64,N_elec-1)) # Gradient of electrolyte concentration
    brug_interface_vector = DiffCache(zeros(Float64,N_elec-1)) # brugeman porosity at the interface
    D_e_node = DiffCache(zeros(Float64,N_elec)) # Electrolyte diffusion coefficient
    D_e_interface = DiffCache(zeros(Float64,N_elec-1)) # Electrolyte diffusion coefficient at the interface
    arg_grad = DiffCache(zeros(Float64,N_elec-1)) # Argument of the gradient term
    grad_ce_flux = DiffCache(zeros(Float64,N_elec)) # Gradient of the electrolyte flux


end 