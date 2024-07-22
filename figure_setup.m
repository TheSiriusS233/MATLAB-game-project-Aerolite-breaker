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
