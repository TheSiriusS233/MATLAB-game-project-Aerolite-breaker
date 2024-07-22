%% function to initialize communication with COM port
% input: com_port = COM port used to communicate with Arduino
% output: arduinoObj = object managing the Serial communication to Arduino
function arduinoObj=init_com(com_port)
    % create a serial client for communication with the Arduino UNO
    arduinoObj=serialport(com_port,115200)
    pause(2); % THIS 2s PAUSE IS CRITICAL FOR MATLAB TO SETUP SERIAL PORT

    % setup I/O terminator to be "Carriage Return" and "Linefeed"
    configureTerminator(arduinoObj,"CR/LF");
    flush(arduinoObj); % discard all data currently in the serial stream
end

