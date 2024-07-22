% Game Project 
% Written by Thursday Group 11
% This code communicates with Arduino to receive joystick data to control
% the initial velocity and angle of the cannonball, after releasing the
% cannonball, Classic RK4 is used to calculate the route of the cannonball, 
% if the cannonball hit the target aerolite, the related score will be added
% and proceed to the next cannonball; else it will directly continue to 
% the next cannonball. The completion of the game is decided by hitting all 
% 5 aerolites using a total of 10 cannonballs or running out all 10 cannonballs.


clear all
close all
clc

com_port='COM3'; % this COM port# should agree with the Arduino IDE 
arduinoObj=init_com(com_port); % initialize communication with Arduino

% define varibles
global mu; % define mu = m_star * G as a global variable used in ODE. 
mu=9*3.986E14; % m^3/s^2, gravitational parameter for the selected star.
r_star=5.024E7; % m, star's radius

% initialize the position of the 5 aerolites
aerolites_x = randi([2E8 5.8E8],1,5); % get random numbers for 5 aerolites' x axis
aerolites_y = randi([-4.3E8 4.3E8],1,5); % get random numbers for 5 aerolites' y axis
aerolite_1 = [aerolites_x(1), aerolites_y(1)]; % giving the first aerolite its position
aerolite_2 = [aerolites_x(2), aerolites_y(2)]; % giving the second aerolite its position
aerolite_3 = [aerolites_x(3), aerolites_y(3)]; % giving the third aerolite its position
aerolite_4 = [aerolites_x(4), aerolites_y(4)]; % giving the fourth aerolite its position
aerolite_5 = [aerolites_x(5), aerolites_y(5)]; % giving the fifth aerolite its position

%% pre-open the interface (figure-window) before communication
[fig,ball,pstar,alphachannel,Background,alphachannel2,...
    aero1,alphachannel3,aero2,alphachannel4,aero3,alphachannel5,aero4,alphachannel6,aero5,alphachannel7,Misaka,alphachannel8] = figure_setup();  %set up the figure window

H_background=image(Background,'Xdata',[-9.1E8 7.9E8],'Ydata',... % define the scale of the background.
                        [-7E8 7E8]);

H_star=image(pstar,'Xdata',[-r_star,r_star],'Ydata',... % define the scale of the star.
                         [-r_star,r_star],'AlphaData',alphachannel);

H_aero1=image(aero1,'Xdata',[aerolite_1(1)-0.5*r_star,aerolite_1(1)+0.5*r_star],'Ydata',... % define the scale of the smallest aerolite.
                         [aerolite_1(2)-0.5*r_star,aerolite_1(2)+0.5*r_star],'AlphaData',alphachannel3);

H_aero2=image(aero2,'Xdata',[aerolite_2(1)-0.625*r_star,aerolite_2(1)+0.625*r_star],'Ydata',... % define the scale of the second smallest aerolite.
                         [aerolite_2(2)-0.625*r_star,aerolite_2(2)+0.625*r_star],'AlphaData',alphachannel4);

H_aero3=image(aero3,'Xdata',[aerolite_3(1)-0.75*r_star,aerolite_3(1)+0.75*r_star],'Ydata',... % define the scale of the third smallest aerolite.
                         [aerolite_3(2)-0.75*r_star,aerolite_3(2)+0.75*r_star],'AlphaData',alphachannel5);

H_aero4=image(aero4,'Xdata',[aerolite_4(1)-0.875*r_star,aerolite_4(1)+0.875*r_star],'Ydata',... % define the scale of the second largest aerolite.
                         [aerolite_4(2)-0.875*r_star,aerolite_4(2)+0.875*r_star],'AlphaData',alphachannel6);

H_aero5=image(aero5,'Xdata',[aerolite_5(1)-r_star,aerolite_5(1)+r_star],'Ydata',... % define the scale of the largest aerolite.
                         [aerolite_5(2)-r_star,aerolite_5(2)+r_star],'AlphaData',alphachannel7);


% print the magnitude and angle of velocity onto the graph
mag_num = 0; % define the initial magnitude as zero
theta = 0; % define the initial angle as zero
score = 0; % define the initial score as zero
ang1 = 0; % store the initial previous rotated angle as zero
frequency = 0; % define the initial frequency as zero

restnum_cannon = 10; % define the initial cannonball number as ten
txt1 = text(-0.07, -0.05, sprintf('Magnitude: %f', mag_num), 'Units', 'normalized'); % display the initial magnitude on the screen
txt2 = text(-0.07, -0.07, sprintf('angle(degree): %f', theta), 'Units', 'normalized'); % display the angle on the screen
txt3 = text(0.65, -0.05, sprintf('score: %.1f', score), 'Units', 'normalized'); % display the score on the screen
txt4 = text(0.65, -0.07, sprintf('rest number of cannon: %d', restnum_cannon), 'Units', 'normalized'); % display the number of cannonball left on the screen
txt7 = text(-0.07, -0.09, sprintf('FPS: %.1f', frequency), 'Units', 'normalized'); % display the fps of the game on the screen
%% initial condition
Matlab_Counter=1;    % initialize Matlab counter for comm synchronization

x_num(1) = 503; % initialize the start of x axis of the joystick as 503, this number is decide by monitoring the arduino monitor
y_num(1) = 503; % initialize the start of y axis of the joystick as 503

Simu_T = 1.75E6;  % simulation time (s)
n=3000;       % times of simulation

% define the gain of initial velocity (ini_V = gain * data from analog pin)
gain_x = 4.45;
gain_y = 4.45;

%% looping area
% The following text gives instructions on how to play this game
txt6 = text(0, 1.02, sprintf('Game Instruction: Adjust the launch Magnitude and the Angle by using joystick and press the button to Launch'),'Units', 'normalized', 'Color', 'Red', 'FontSize',14);   

for num_of_cannon = 1:10  % 1 connon balls per loop, a total of 10 loops
    x=zeros(4,n+1); % preallocate RK4 Method solution array for speed
    t=zeros(n+1,1); % preallocate time vector for speed
    n=2204;       % refresh the number of n

    while 1 % calculating the movement of the joystick, continue to loop until the button is pressed
        Misaka = imrotate(Misaka, -ang1,'crop'); % reverse rotate the amount that the railgun rotate in the previous loop
        writeline(arduinoObj,int2str(Matlab_Counter)); % send Matlab counter to Arduino (used for comm synchronization)
        ang_image = 0; % preset the angle of the railgun to be 0

        % wait for and receive counter and joystick data from Arduino
        data=data_from_Arduino(arduinoObj);
        tmp=split(data,',');    % data from Arduino is comma-delimited
        num=str2double(tmp);    % convert strings to numeric double type array
        Arduino_Counter=num(1); % store Arduino counter data (synchronization)
    
        x_num(Matlab_Counter+1) = num(2); % store the x axis of the joystick
        y_num(Matlab_Counter+1) = num(3); % store the y axis of the joystick

        if num(2) == 504         % Once detect num(2) - 504 = 0, add num(2) by 0.1
            num(2) = num(2)+0.1; % to avoid (num(2)-504)*gain_x) becomes 0. (num/0 is invaid)
        end

        mag_num = sqrt(((num(2)-504)*gain_x)^2 + ((num(3)-504)*gain_y)^2); % calculate the magnitude of the velocity
        theta = atan(((num(3)-504)*gain_y)/((num(2)-504)*gain_x))/pi*180; % calculate the angle of the velocity
        if 90>theta&&theta>0&&num(2)-504<0||0>theta&&theta>-90&&num(2)-504<0 % apply the angle to the negative side of the x axis
            theta=theta+180; % renew the angle if the velocity of the x axis is negative
        end
        ang_image = ang_image + theta; % apply the velocity angle to the angle of the railgun
        Misaka = imrotate(Misaka, -ang_image,'crop'); % rotate the railgun according to the angle calculated
        ang1 = -ang_image; % store the angle rotated to ang1
        H_railgun=image(Misaka,'Xdata',[-33*0.5*r_star,-27*0.5*r_star],'Ydata',... % define the scale of the ball.
                         [-3*0.5*r_star,3*0.5*r_star],'AlphaData',alphachannel8);
        set(H_railgun,'Clipping','off'); % prevent image from cropping at edge of axes
        
        % check synchronization problem
        if Arduino_Counter~=Matlab_Counter
            % if counters don't match then comm synchronization problem
            disp('Communication synchronization problem, stopping program.') 
            return  % stop execution of program
        end
        
        % compare the last result to the current result (release the cannonball).
        if num(4) ~= 1 % if the button is pressed
            break % jump out of the while loop
        end
        Matlab_Counter = Matlab_Counter + 1; % update Matlab Counter values to check for synchronization 
        % update text object
        set(txt1, 'String', sprintf('Magnitude: %f', mag_num)); % update the magnitude of the velocity
        set(txt2, 'String', sprintf('angle(degree): %f', theta)); % update the angle of the velocity
        
        pause(0.08);
        delete(H_railgun);
    end   

    Matlab_Counter = Matlab_Counter + 1; % update Matlab Counter value for the final time to extract the final value of the velocity
    t0=0; tn=Simu_T; % time interval from [t0 to tn] (s)
     % n time steps to use for solution
    h=(tn-t0)/n;% s, step size
     % set up matrix
    x0=[-15*r_star;             % x_position
        0;                      % y_position
        (gain_x*(x_num(Matlab_Counter)-510));  % x_velocity
        (gain_y*(y_num(Matlab_Counter)-510))]; % y_velocity
    restnum_cannon = restnum_cannon - 1; % decrease the number of cannonball by 1
    set(txt4, 'String', sprintf('rest number of cannon: %.1d', restnum_cannon)); % update the number of cannonball left

    % graph plotting loop 
    for i=1:n % main iteration loop
        tic;  % start to measure time per loop
        %RK4 calculation to get position and velocity of each n.
        w1=1/6; w2=1/3; w3=1/3; w4=1/6; c2=1/2; c3=1/2; c4=1;
        a21=1/2; a31=0; a32=1/2; a41=0; a42=0; a43=1;
        t(1)=t0; % set initial condition for time
        x(:,1)=x0; % set initial condition for solution
        % increment i from 1 to n, note n+1 points in solution
        x(:,i+1)=RK4(t(i),x(:,i),h,w1,w2,w3,w4,c2,c3,c3,a21,...
        a31,a32,a41,a42,a43); % Classic RK4 Method
        t(i+1)=t(i)+h; % increment time by the time step h (s)
    
        % plot initial ball location in figure window
        
        [a,b,~]=size(ball); % determine size of ball image, ~ means ignore the third term (color)
        scale=0.87E5;         % set scaling of ball image in figure
     
        H=image(ball,'Xdata',[x(1,i)-scale*a/2,x(1,i)+scale*a/2],'Ydata',... % define the scale of the ball.
                             [x(2,i)-scale*b/2,x(2,i)+scale*b/2],'AlphaData',alphachannel2); % row 1 is position, row 2 is velocity, column is the corresponing time.
        set(H,'Clipping','off'); % prevent image from cropping at edge of axes
    
         % delect the cannonball and begin another turn if it is out of boundary.
         % x-axis
         if x(1,i) > 8.9E8 || x(1,i) < -8.9E8   % range of x
            delete(H) % delete the image of the cannonball
            break; % jump out of the loop
         end
         % y-axis
         if x(2,i) > 6.65E8 || x(2,i) < -6.8E8    % range of y
            delete(H) % delete the image of the cannonball
            break; % jump out of the loop
         end
         % delete the cannonball when impact the star.
         if sqrt( (x(1,i))^2 + (x(2,i))^2 ) < r_star*1.09  % radius of the star+cannonball
            delete(H) % delete the image of the cannonball
            break; % jump out of the loop
         end

         % detect whether the cannonball hit the aerolite
         if x(1,i) > aerolite_1(1)-0.5*r_star && x(1,i) < aerolite_1(1)+0.5*r_star &&...
          x(2,i) > aerolite_1(2)-0.5*r_star  && x(2,i) < aerolite_1(2)+0.5*r_star % aerolite_1
            
            aerolite_1(1) =  15E8; % use the number out of boundary to 
            aerolite_1(2) = -10E8; % disable this if once the aerosite is hitted
            score = score + 3; % update the score
            delete(H_aero1) % delete the image of the aerolite
            delete(H) % delete the image of the cannonball
            break % jump out of the loop
         end

          if x(1,i) > aerolite_2(1)-0.625*r_star*1.12 && x(1,i) < aerolite_2(1)+0.625*r_star*1.12 &&...
          x(2,i) > aerolite_2(2)-0.625*r_star*1.12  && x(2,i) < aerolite_2(2)+0.625*r_star*1.12 % aerolite_2

            aerolite_2(1) =  15E8; % use the number out of boundary to 
            aerolite_2(2) = -10E8; % disable this if once the aerosite is hitted
            score = score + 2.5; % update the score
            delete(H_aero2) % delete the image of the aerolite
            delete(H) % delete the image of the cannonball
            break % jump out of the loop
          end

           if x(1,i) > aerolite_3(1)-0.75*r_star*1.12 && x(1,i) < aerolite_3(1)+0.75*r_star*1.12 &&...
          x(2,i) > aerolite_3(2)-0.75*r_star*1.12  && x(2,i) < aerolite_3(2)+0.75*r_star*1.12 % aerolite_1

            aerolite_3(1) =  15E8; % use the number out of boundary to 
            aerolite_3(2) = -10E8; % disable this if once the aerosite is hitted
            score = score + 2; % update the score
            delete(H_aero3) % delete the image of the aerolite
            delete(H) % delete the image of the cannonball
            break % jump out of the loop
           end

            if x(1,i) > aerolite_4(1)-0.875*r_star*1.12 && x(1,i) < aerolite_4(1)+0.875*r_star*1.12 &&...
          x(2,i) > aerolite_4(2)-0.875*r_star*1.12  && x(2,i) < aerolite_4(2)+0.875*r_star*1.12 % aerolite_1

            aerolite_4(1) =  15E8; % use the number out of boundary to 
            aerolite_4(2) = -10E8; % disable this if once the aerosite is hitted
            score = score + 1.5; % update the score
            delete(H_aero4) % delete the image of the aerolite
            delete(H) % delete the image of the cannonball
            break % jump out of the loop
           end

            if x(1,i) > aerolite_5(1)-1.12*r_star && x(1,i) < aerolite_5(1)+1.12*r_star &&...
          x(2,i) > aerolite_5(2)-1.12*r_star  && x(2,i) < aerolite_5(2)+1.12*r_star % aerolite_1

            aerolite_5(1) =  15E8; % use the number out of boundary to 
            aerolite_5(2) = -10E8; % disable this if once the aerosite is hitted
            score = score + 1; % update the score
            delete(H_aero5) % delete the image of the aerolite
            delete(H) % delete the image of the cannonball
            break % jump out of the loop
           end
         
        pause(0.001); % wait for ball to display,faster than drawnow() function
        delete(H) % delete the previous image plot.
    
        frequency = 1/toc; % calculate the frequency
        set(txt7, 'String', sprintf('FPS: %.1f', frequency)); % display the fps of the game of the screen
    end

    set(txt3, 'String', sprintf('score: %.1f', score)); % display the score of the game of the screen

    if score == 10 % determine if full score is reached
        score = score + restnum_cannon*2.5; % add bonus points, 2.5 points per cannonball left
        set(txt3, 'String', sprintf('score: %.1f', score)); % update the score on the screen
        break % jump out of the loop
    end

    if num_of_cannon == 1 % delete the game instruction after the first cannonball run
        delete(txt6); % delete the instruction
    end
end

txt5 = text(-465000000, 0, sprintf('Your Final Score: %.1f', score), 'Color', 'Yellow', 'FontSize',50); % display the final score on the center of the screen

%% function to setup figure window
% inputs: none
% outputs: f = handle to figure window (provides a reference to the figure)
%       ball = image of black ball (flipped so it is oriented correctly)

function [f,ball,pstar,alphachannel,Background,alphachannel2,aero1,alphachannel3,aero2,alphachannel4,aero3,alphachannel5...
    ,aero4,alphachannel6,aero5,alphachannel7,Misaka,alphachannel8]=figure_setup()
    % adjust figure to desired size&position, then use f.Position to get #s
    f=figure('position',[10     42    1600    1000]);   % B316, figure on LHS of screen
    set(gcf,'Toolbar','none','Menu','none');    % remove toolbar and menu
    set(gcf,'color','w');     % make figure background white
    set(gca,'visible','off'); % remove axis labels
    hold on        % allow superposition of plots
    xlim([-9.1E8 7.9E8]) % set-up x-limit
    axis('equal')  % set aspect ratio so equal tick mark increments on each 
                   % axis are equal in size
    axis('manual') % freeze all axis limits for subsequent plots so they 
                   % do not automatically adjust on the fly
    % read image of the star to track
    [pstar,~,alphachannel] = imread('star.png'); 
    pstar = flipud(pstar); % need to flip image so it is oriented correctly
    alphachannel = flipud(alphachannel); % remove the part outside the edges of the image

    % read image of the background to track
    Background = imread('background.png'); % need to flip image so it is oriented correctly
    Background = flipud(Background); % remove the part outside the edges of the image

    % read image of the cannonball to track
    [ball,~,alphachannel2] = imread('cannon.png'); 
    ball=flipud(ball); % need to flip image so it is oriented correctly
    alphachannel2 = flipud(alphachannel2); % remove the part outside the edges of the image
    
    % read image of the aerolite to track
    [aero1,~,alphachannel3] = imread('aerolite.png'); 
    aero1=flipud(aero1); % need to flip image so it is oriented correctly
    alphachannel3 = flipud(alphachannel3); % remove the part outside the edges of the image

    % read image of the aerolite to track
    [aero2,~,alphachannel4] = imread('aerolite.png'); 
    aero2=flipud(aero2); % need to flip image so it is oriented correctly
    alphachannel4 = flipud(alphachannel4); % remove the part outside the edges of the image

    % read image of the aerolite to track   
    [aero3,~,alphachannel5] = imread('aerolite.png'); 
    aero3=flipud(aero3); % need to flip image so it is oriented correctly
    alphachannel5 = flipud(alphachannel5); % remove the part outside the edges of the image
    
    % read image of the aerolite to track
    [aero4,~,alphachannel6] = imread('aerolite.png'); 
    aero4=flipud(aero4); % need to flip image so it is oriented correctly
    alphachannel6 = flipud(alphachannel6); % remove the part outside the edges of the image

    % read image of the aerolite to track
    [aero5,~,alphachannel7] = imread('aerolite.png'); 
    aero5=flipud(aero5); % need to flip image so it is oriented correctly
    alphachannel7 = flipud(alphachannel7); % remove the part outside the edges of the image

    % read image of the cannon to track
    [Misaka,~,alphachannel8] = imread('railgun.png'); 
    Misaka=flipud(Misaka); % need to flip image so it is oriented correctly
    alphachannel8 = flipud(alphachannel8); % remove the part outside the edges of the image
end


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

%% function containing the differential equation dx/dt=f(t,x) to solve
% inputs:     t = time (s)
%             x = state vector (4x1) consists of mass position and velocity
function dxdt=f(t,x)
global mu;
r=norm(x(1:2),2); % magnitude of position vector
dxdt=[zeros(2,2) eye(2);
-(mu/r^3)*eye(2) zeros(2,2)]*x; % dxdt is a 4x4 vector
end

%% function to implement Runge-Kutta (RK) 4th order method
% inputs: ti = current time (s)
%         xi = previous state vector estimate (4x1): position and velocity
%          h = sample time (s)
%  w1,w2,w3,w4,c2,c3,c4,a21,a31,a32,a41,a42,a43 = "Classic" RK4 parameters
% outputs: x_new = new state vector estimate (4x1): position and velocity
function x_new=RK4(ti,xi,h,w1,w2,w3,w4,c2,c3,c4,...
a21,a31,a32,a41,a42,a43)
k1=h*f(ti,xi); % output of f(t,x) is 4x1, so k1 is 4x1
k2=h*f(ti+c2*h,xi+a21*k1);
k3=h*f(ti+c3*h,xi+a31*k1+a32*k2);
k4=h*f(ti+c4*h,xi+a41*k1+a42*k2+a43*k3);
x_new=xi+w1*k1+w2*k2+w3*k3+w4*k4;
end
