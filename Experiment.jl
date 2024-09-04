module Experiment_func

    using MAT

    function get_capacity_data(temperature::Int, soc::Int)
        data = matread("RPT_analysis_data.mat")
        specific_data = get(get(data, "Temperature_$temperature", Dict()), "SOC_$soc", nothing)
        isnothing(specific_data) && error("Data not found for Temperature $temperature and SOC $soc")
        return specific_data["Days"], specific_data["Capacity_Mean"], specific_data["Capacity_Std"]
    end

    function get_LAM_data(temp::Int, soc::Int)
        data = matread("RPTx_analysis_data.mat")
        specific_data = get(get(data, "Temperature_$temp", Dict()), "SOC_$soc", nothing)
        isnothing(specific_data) && error("No data for Temperature $temp and SOC $soc")
        return specific_data["LAM_days"], specific_data["LAM_mean"], specific_data["LAM_std"]
    end

    function Storage_step(storage_SOC,Storage_duration_btn_RPT_days)

        Time_to_storage_SOC = (4.8/1.6)*3600.0*(1.0-storage_SOC)
        Storage_duration_btn_RPT  = Storage_duration_btn_RPT_days*24*3600.0 
        ageing_exp = [(false, x -> 0.0, "Time", 3600.0, 5000.0, nothing),
        (false, x -> -1.6, "Voltage", 4.2, 3600 * 10, nothing),
        (true, x -> -1.6, "Abs_Current", 0.25, 50000.0, 4.2),
        (false, x -> 0.0, "Time", 1800.0, 50000.0, nothing),
        (false, x -> 1.6, "Time", Time_to_storage_SOC, Time_to_storage_SOC * 2.0, nothing),
        (false, x -> 0.0, "Time", Storage_duration_btn_RPT, Storage_duration_btn_RPT * 2.0, "Storage_step"),
        (false, x -> 1.6, "Voltage", 2.5, 3600.0 * 5.0, nothing)]

        return ageing_exp

    end 

    function Calendar_ageing_exp_from_dates(dates,storage_SOC)


        Discharge_to_zero_SOC_step = (false,x->1.6,"Voltage", 2.5,3600*1.1,nothing,600.0)

        RPT_cycle = [(false,x->0.0,"Time", 3600.0,5000.0,nothing,600.0),
        (false,x->-1.6,"Voltage", 4.2,3600*10,nothing,300.0),
        (true,x->-1.6,"Abs_Current", 0.25,50000,4.2,180.0),
        (false,x->0.0,"Time", 1800,50000,nothing,400.0),
        (false,x->1.6,"Voltage", 2.5,3600*5.0,"RPT-Capacity",300.0)]

        ageing_exp = [Discharge_to_zero_SOC_step;RPT_cycle]

        for i in 1:lastindex(dates)-1
            Storage_duration_btn_RPT_days = dates[i+1] - dates[i]
            ageing_exp = [ageing_exp;Storage_step(storage_SOC,Storage_duration_btn_RPT_days);RPT_cycle]
        end
        return ageing_exp

    end

    function Storage_step(storage_SOC,Storage_duration_btn_RPT_days,Capacity)

        Time_to_storage_SOC = (Capacity/1.6)*3600.0*(1.0-storage_SOC)
        Storage_duration_btn_RPT  = Storage_duration_btn_RPT_days*24*3600.0 
        if storage_SOC ==1.0

            ageing_exp = [(false, x -> 0.0, "Time", 3600.0, 5000.0, nothing,600.0),
            (false, x -> -1.6, "Voltage", 4.2, 3600 * 10, nothing,300.0),
            (true, x -> -1.6, "Abs_Current", 0.25, 50000.0, 4.2,180.0),
            (false, x -> 0.0, "Time", 1800.0, 50000.0, nothing,400.0),
            (false, x -> 0.0, "Time", Storage_duration_btn_RPT, Storage_duration_btn_RPT * 2.0, "Storage_step",Storage_duration_btn_RPT/5.0),
            (false, x -> 1.6, "Voltage", 2.5, 3600.0 * 5.0, nothing,300.0)]

        else
            ageing_exp = [(false, x -> 0.0, "Time", 3600.0, 5000.0, nothing,600.0),
            (false, x -> -1.6, "Voltage", 4.2, 3600 * 10, nothing,300.0),
            (true, x -> -1.6, "Abs_Current", 0.25, 50000.0, 4.2,180.0),
            (false, x -> 0.0, "Time", 1800.0, 50000.0, nothing,400.0),
            (false, x -> 1.6, "Time", Time_to_storage_SOC, Time_to_storage_SOC * 2.0, nothing,Time_to_storage_SOC/5.0),
            (false, x -> 0.0, "Time", Storage_duration_btn_RPT, Storage_duration_btn_RPT * 2.0, "Storage_step",Storage_duration_btn_RPT/5.0),
            (false, x -> 1.6, "Voltage", 2.5, 3600.0 * 5.0, nothing,300.0)]
        end 

        return ageing_exp

    end 

    function Calendar_ageing_exp_from_dates(dates,storage_SOC,Exp_capacity)


        Discharge_to_zero_SOC_step = (false,x->1.6,"Voltage", 2.5,3600*1.1,nothing,600.0)

        RPT_cycle = [(false,x->0.0,"Time", 3600.0,5000.0,nothing,600.0),
        (false,x->-1.6,"Voltage", 4.2,3600*10,nothing,300.0),
        (true,x->-1.6,"Abs_Current", 0.25,50000,4.2,180.0),
        (false,x->0.0,"Time", 1800,50000,nothing,400.0),
        (false,x->1.6,"Voltage", 2.5,3600*5.0,"RPT-Capacity",300.0)]

        ageing_exp = [Discharge_to_zero_SOC_step;RPT_cycle]

        for i in 1:lastindex(dates)-1
            Capacity = Exp_capacity[i]
            Storage_duration_btn_RPT_days = dates[i+1] - dates[i]
            ageing_exp = [ageing_exp;Storage_step(storage_SOC,Storage_duration_btn_RPT_days,Capacity);RPT_cycle]
        end
        return ageing_exp

    end


end 