%% function to receive data from Arduino when it is available
% input: arduinoObj = object managing the Serial communication to Arduino
% output:      data = comma-delimited data from Arduino
function data=data_from_Arduino(arduinoObj)
    % wait for Arduino to have data available for Matlab to read
    while true
        if arduinoObj.NumBytesAvailable>0 % Arduino has data for Matlab
            break; % exit while true loop
        end
    end
    % receive comma-delimited data from Arduino
    data=readline(arduinoObj);
end

