function strRes = t_res_geodata(time,timeUnits)

if numel(time) == 1
    strRes = 'unknown';
else
    dTime = mean(diff(time));

    if regexpbl(timeUnits,'day')
        if dTime > 20 && dTime < 35
            strRes = 'month';
        elseif dTime > 0.9 && dTime < 1.1
            strRes = 'day';
        else
            error('t_res_geodata:timeStep',['The average time step is ' num2str(dTime) ' does not correspond to a specific time step.']); 
        end
    elseif regexpbl(timeUnits,'hour')
        if dTime > 22 && dTime < 26
            strRes = 'day';
        elseif dTime < 22
            strRes = 'hour';
        elseif dTime > 600 && dTime < 750
            strRes = 'month';
        else
            error('t_res_geodata:timeStep',['The average time step is ' num2str(dTime) ' does not correspond to a specific time step.']); 
        end
    else
       error('t_res_geodata:unknownUnits',['The units ' char(39) timeUnits char(39) ' are not recognized.']); 
    end
end

    
    