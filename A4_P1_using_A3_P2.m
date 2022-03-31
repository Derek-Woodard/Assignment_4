%%
% Derek Woodard, 100827083

%% Part 2

clear all
close all
set(0,'DefaultFigureWindowStyle','docked');

m0 = 9.10938356e-31;            % rest mass of an electron
mn = 0.26*m0;                   % effective mass of electrons
T = 300;                        % temperature
k = 1.38064852e-23;             % Boltzmann's Constant
Tmn = 0.2e-12;                  % mean time between collisions
q = 1.60217662e-19;             % The charge of an electron
vth = sqrt((2*k*T)/mn);         % calculate thermal velocity for 1.1 with two degrees of freedom
mfp = vth * Tmn;                % calculate the mean free path for 1.2


%% The nominal size of the region is 200nm x 100nm
H = 100e-9;
W = 200e-9;

pop_tot = 1500;                 % use this variable to determine how many electrons are modeled
pop_vis = 5;                    % use this variable to change how many electrons are plotted
ts = H/vth/100;                 % The time step is set to make each plot point have an electron move 1/100th of the height of the area
iter = 3000;                    % how many plot points total
Vx = 0.8;                       % voltage across the x direction is set to 0.8V
Vy = 0;                         % voltage across the y dimension not provided, so assume 0
m0 = 9.10938215e-31;            % mass of an electron
density = 1e19;                 % density of electrons
ps = 1 - exp(-((ts)/(Tmn)));    % The eletrons have a probability of scattering

animate = 1;                    % set to 1 to see each iteration
% animate = 0;                    % set to 0 to only see final product

initvels = zeros(pop_tot,1);
% The initially set velocities of all atoms are stored here to be tracked on
% a histogram

mfps = zeros(pop_tot,1);
meanfp = 0;
% To measure the mean free path, an array is used to track how far each
% electron has travelled before scattering

tbc = zeros(pop_tot,1);
meantbc = 0;
% The measured time between collisions will be tracked through an array
    
% Set applied voltage to 0.8
Vo = 0.8;

% For scale simplicity, the nano-factor is removed - This will be
% re-applied after the acceleration is found
scale = 1e-9;
H = round(H/scale);
W = round(W/scale);
ts = ts/1e-9;

% Set the spacing for the mesh and the number of nodes total based on
% region size
dx = 1;
dy = 1;
nx = round(W/dx);
ny = round(H/dy);

% There are two different sigma values based on the bottleneck
sig1 = 1;
sig2 = 10^-2;
sigs = zeros(nx,ny);

% We need to set the 'bottleneck' by creating two boxes in the region
boxes = [nx*2/5 nx*3/5 ny*2/5 ny*3/5];

% A circle is formed by filling the following: [center_xpos center_ypos radius]
% circ = scale .* [150 50 10];
% theta = linspace(0,2*pi);
% xc = circ(3)*cos(theta) + circ(1);
% yc = circ(3)*sin(theta) + circ(2);
% The circle is only used in the trajectory portion, not the potential
% calculation

% We need to set the sigma values for the full region, giving different
% sigma values to the boxes
for i = 1:nx
    for j = 1:ny
        if i > boxes(1) && i < boxes(2) && (j < boxes(3) || j > boxes(4))
            sigs(i,j) = sig2;          
        else
            sigs(i,j) = sig1;
        end
    end
end

%%
% The G and F matrices are formed using the number of points calculated above
G = zeros(nx*ny);
F = zeros(1,nx*ny);

for i=1:nx
    for j=1:ny
        % set up the mapping variables
        n = j+(i-1)*ny;
        nxp = j+i*ny;
        nxm = j+(i-2)*ny;
        nyp = j+1+(i-1)*ny;
        nym = j-1+(i-1)*ny;
        
        if i == 1
            G(n,:) = 0;
            G(n,n) = 1;
            F(n) = Vo;
        elseif i == nx
%             G(n,:) = 0;
            G(n,n) = 1;
%             F(n) = 0;
        elseif j == 1
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp));
        elseif j == ny
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nym));
        else
            G(n,nxp) = (sigs(i+1,j) + sigs(i,j))/2;
            G(n,nxm) = (sigs(i-1,j) + sigs(i,j))/2;
            G(n,nyp) = (sigs(i,j+1) + sigs(i,j))/2;
            G(n,nym) = (sigs(i,j-1) + sigs(i,j))/2;
            G(n,n) = -(G(n,nxp) + G(n,nxm) + G(n,nyp) + G(n,nym));
        end
    end
end

V = G\F';
S = zeros(ny,nx,1);

for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        S(j,i) = V(n);
    end
end

figure(1)
mesh(S)
colormap(jet)
xlabel('x')
ylabel('y')
title('Bottleneck Potential')

[Exs Eys] = gradient(S);
Eys = -Eys;
Exs = -Exs;

figure(2)
quiver(Exs, Eys)
ylabel('y');
xlabel('x');
title('Electric Field with Bottleneck')

Exs = Exs';
Eys = Eys';
Fxs = zeros(W,H);
Fys = zeros(W,H);
Axs = zeros(W,H);
Ays = zeros(W,H);

% Since the electric field uses numbered blocks instead of measurements, we
% can return the scale to normal and still get the same results
H = H*scale;
W = W*scale;
ts = ts*scale;
boxes = scale .* boxes;

%%
% In order to track the electrons, the following arrays will be used
% The arrays will store the positions, velocities, and temperatures of the
% elctrons
state = zeros(pop_tot, 4);
traj = zeros(iter, pop_vis*2);
temps = zeros(iter,1);
avgtemps = zeros(iter,1);
J = zeros(iter,2);

% %% 
% % Start by generating the initial eletrons with a constant speed of vth
% for i = 1:pop_tot
%     angle = rand*2*pi;
%     vel = randn*vth;
%     initvels(i) = vel;
%     state(i,:) = [W*rand H*rand vel*cos(angle) vel*sin(angle)];
%         while(state(i,1) > boxes(1) && state(i,1) < boxes(2) && (state(i,2) < boxes(3) || state(i,2) > boxes(4)))
%             state(i,1:2) = [W*rand H*rand];
%         end
% end
% This was removed because the electons are being injected in from the left
% side.




%%
% Using the calculated electric field, we determine the force, then
% acceleration at every 1x1 location
Fxs(:,:) = q.* Exs(:,:);
Fys(:,:) = q.* Eys(:,:);
Axs(:,:) = Fxs(:,:)./m0;
Ays(:,:) = Fys(:,:)./m0;
% Given that these values are derived for incredibly small areas, the
% resulting accelerations are not noticable when plotting the trajectories.
% For this reason, these values are not used in the plotting, however, the
% code that uses them is commented out to show how they can be used.
% Ex = Vx/W;
% Ey = Vy/H;
% Fx = q * Ex;
% Fy = q * Ey;
% Ax = Fx/m0;
% Ay = Fy/m0;

xloc = zeros(pop_tot,1);
yloc = zeros(pop_tot,1);

%%
% Now we iterate to update the electron positions and plot the trajectories
for i = 1:iter
    % Inject the electrons from the left side
    if i <= pop_tot
        angle = rand*2*pi;
        vel = randn*vth;
        initvels(i) = vel;
        state(i,:) = [0 ((H*0.6 - H*0.4)*rand)+(H*0.4) abs(vel*cos(angle)) vel*sin(angle)];
    end
  
    state(:,1:2) = state(:,1:2) + ts.*state(:,3:4);
    % update the state of each electron using the speed
    
    % We need to ensure the boundary conditions are checked
    bcx = state(:,1) > W;
%     state(bcx,1) = state(bcx,1) - W;
    % Ensure the electron is transported from the right side to the left if
    % travelling past the right X limit
    state(bcx,1) = 0;
    state(bcx,2) = ((H*0.6 - H*0.4)*rand)+(H*0.4);
    state(bcx,3) = abs(randn*vth*cos(rand*2*pi));
    state(bcx,4) = randn*vth*sin(rand*2*pi);
    % The above checks to see if an electron has escaped to the left or the
    % right. If it has, it is reset to the initial injection location and
    % given a new speed. This ensures that there is no loss of electrons
    % while maintaining the distribution.
    
    bcx = state(:,1) < 0;
%     state(bcx,1) = state(bcx,1) + W;
    % Ensure the electron is transported from the left side to the right if
    % travelling past the left X limit
    state(bcx,1) = 0;
    state(bcx,2) = ((H*0.6 - H*0.4)*rand)+(H*0.4);
    state(bcx,3) = abs(randn*vth*cos(rand*2*pi));
    state(bcx,4) = randn*vth*sin(rand*2*pi);
    
    bcy = state(:,2) > H;
    state(bcy,2) = 2*H - state(bcy,2);
    state(bcy,4) = -state(bcy,4);
    % Ensure the electrons bounce back down if going above high the Y limit
    
    bcy = state(:,2) < 0;
    state(bcy,2) = -state(bcy,2);
    state(bcy,4) =- state(bcy,4);
    % Ensure the electrons bounce back up if going below the low Y limit
    
    % Check for collisions with the boxes
    for b = 1:pop_tot
        while(state(b,1) > boxes(1) && state(b,1) < boxes(2) && (state(b,2) < boxes(3) || state(b,2) > boxes(4)))
            if state(b,3) > 0
                xd = state(b,1) - boxes(1);
                newx = boxes(1);
            else
                xd = boxes(2) - state(b,1);
                newx = boxes(2);
            end
            
            if(state(b,4) > 0)
                yd = state(b,2) - boxes(4);
                newy = boxes(4);
            else
                yd = boxes(3) - state(b,2);
                newy = boxes(3);
            end
            
            if(xd < yd)
                state(b,1) = newx;
                state(b,3) = -state(b,3);
            else
                state(b,2) = newy;
                state(b,4) = -state(b,4);
            end
        end
    end
    
    % Check for collisions with the circle
%     for b = 1:pop_tot
%         if (state(b,1)-circ(1))^2 + (state(b,2)-circ(2))^2 < (circ(3))^2
%             nor = atan2(state(b,2)-circ(2),state(b,1)-circ(1));
% 
%             vpart = [state(b,3) state(b,4)];
%             vn = norm(vpart);
% 
%             circlecenterv = [state(b,1)-circ(1) state(b,2)-circ(2)];
%             ccn = norm(circlecenterv);
% 
%             incangle = acos(dot(-vpart,circlecenterv)/(vn*ccn));
% 
%             new_angle = nor - incangle;
% 
%             v = sqrt((state(b,3))^2 + (state(b,4))^2);
% 
%             vx = cos(new_angle)*v;
%             vy = sin(new_angle)*v;
% 
%             state(b,3) = vx;
%             state(b,4) = vy;
%             state(b,1) = state(j,1) + ts*state(j,3);
%             state(b,2) = state(j,2) + ts*state(j,4);
%         end
%     end
    
        % This code uses the acceleration matrices that do not provide
        % enough to visibly affect the plot
    % Determine where in the area the electons are
    xloc(:) = round(state(:,1));
    xloc(xloc==0) = 1;
    
    yloc(:) = round(state(:,2));
    yloc(yloc==0)=1;
    
    for b = 1:pop_tot:1
        state(b,3) = state(b,3) + ts*Axs(xloc(b),yloc(b));
        state(b,4) = state(b,4) + ts*Ays(xloc(b),yloc(b));
    end
%     % The speed of each electron is updated using the acceleration from the
%     % electric field applied
%     state(:,3) = state(:,3) + ts*Ax;
%     state(:,4) = state(:,4) + ts*Ay;
    
    temps(i) = (sum(state(:,3).^2) + sum(state(:,4).^2)) * mn/k/2/pop_tot;
    avgtemp = sum(temps(1:i,1))/i;
    avgtemps(i) = avgtemp;
    % Track the temperature using the velocity of every electron at every
    % point
    % Track the temperature using the velocity of every electron at every
    % point
    
    for j = 1:pop_vis
        traj(i, (2*j):(2*j+1)) = state(j, 1:2);
    end
    % Track the trajectory of a small number of electrons by storing the x
    % and y position for each iteration
    
    ra = rand(pop_tot, 1) < ps;
    scattersize = size(state(ra,3:4));
    distvels = zeros(scattersize(1),1) + (randn(scattersize(1),1).*vth);
    newvels = ones(scattersize);
    
    for a = 1:scattersize(1)-1
        angle = rand*2*pi;
        vel = randn*vth;
        newvels(a, :) = [vel*cos(angle) vel*sin(angle)];
    end
    
    mfps(:) = mfps(:) + (ts * (sqrt(state(:,3).^2 + state(:,4).^2)));
    mfps(ra) = 0;
    meanfp = sum(mfps)/pop_tot;
    % The measured mean free path is found by calculating the distance
    % between each scatter for all electrons
    
    tbc(:) = tbc(:) + ts;
    tbc(ra) = 0;
    meantbc = sum(tbc)/pop_tot;  
 
    state(ra,3:4) = newvels;
    % Check to see if the electrons will scatter based around the
    % probability calculated
    % If they do scatter, the velocity and direction will be changed to a
    % random value within the normal distribution around vth
    
    % The animation does not need to be updated at every iteration
    % So we update it every 5 iterations to use less memory
    if animate && mod(i,5) == 0
        figure(4)
        pbaspect([2,1,1]);
        subplot(2,1,1);
        hold off;
        plot(state(1:pop_vis,1), state(1:pop_vis,2), 'o');
        title('Electrons an Electric Field Applied (P.2)');
        xlabel('X (nm)');
        ylabel('Y (nm)');
        hold on;
        for j=1:pop_vis
            plot(traj(:,j*2), traj(:,j*2+1), '.');
        end
        
        % Plot the boxes
        rectangle('Position',  [boxes(1) boxes(4) boxes(2)-boxes(1) (scale*ny)-boxes(4)]);
        rectangle('Position',  [boxes(1) 0 boxes(2)-boxes(1) (scale*ny)-boxes(4)]);
        
        % Plot the circle
%         plot(xc,yc);
%         hold on;
        
%         if i > 1
%             subplot(2,1,2);
%             hold off;
%             plot(ts*(0:i-1), temps(1:i));
%             hold on
%             plot(ts*(0:i-1), avgtemps(1:i));
%             axis([0 ts*iter min(temps)*0.98 max(temps)*1.02]);
%             title(sprintf('Temperature (avg = %f)', avgtemp));
%             xlabel('Time (s)');
%             ylabel('Temperature (K)');
%         end
        pause(0.05)
    end
end

%%
% We track the distribution of the speeds using a histogram
% The average value of the initial values should be around vth
% vavg = sqrt(sum(state2(:,3).^2)/pop_tot + sum(state2(:,4).^2)/pop_tot);
% figure(3)
% histogram(initvels)
% title(sprintf('Initial Velocities - avg = %d', vavg));
% xlabel('velocity (m/s)');
% ylabel('electrons');


% Now we need to display the trajectory after the iterations have completed
figure(4);
pbaspect([2,1,1]);
subplot(2,1,1);
title('Electrons an Electric Field Applied (P.2)');
xlabel('X (nm)');
ylabel('Y (nm)');
hold on;
for i=1:pop_vis
    plot(traj(:,i*2), traj(:,i*2+1), '.');
end

% Plot the boxes
rectangle('Position',  [boxes(1) boxes(4) boxes(2)-boxes(1) (scale*ny)-boxes(4)]);
rectangle('Position',  [boxes(1) 0 boxes(2)-boxes(1) (scale*ny)-boxes(4)]);

% Plot the circle
% plot(xc,yc);
% hold on;

% if i > 1
%     subplot(2,1,2);
%     hold off;
%     plot(ts*(0:iter-1), temps(1:iter));
%     hold on
%     plot(ts*(0:iter-1), avgtemps(1:iter));
%     axis([0 ts*iter min(temps)*0.98 max(temps)*1.02]);
%     title(sprintf('Temperature (avg = %f)', avgtemp));
%     xlabel('Time (s)');
%     ylabel('Temperature (K)');
% end

% % An electron density map is created for the bottlenecked simulation
% dmap = [state(:,1), state(:,2)];
% figure(5)
% hist3(dmap, [200 100]);
% surfHandle = get(gca, 'child');
% set(surfHandle, 'FaceColor', 'interp', 'EdgeColor', 'none', 'CdataMode', 'auto');
% view(2)
% title('Electron Density Map');
% xlabel('x (m)');
% ylabel('y (m)');


%%
% Similarily to the density map, a temperature map is created
% This map is created by putting the electron velocities into areas much
% smaller than the full space called bins.
% x_temp_sum = zeros(ceil(W/1e-9), ceil(H/1e-9));
% y_temp_sum = zeros(ceil(W/1e-9), ceil(H/1e-9));
% num_temp = zeros(ceil(W/1e-9), ceil(H/1e-9));

% Start by observing the velocities of all the electrons
% for i = 1:pop_tot
%     % We now need to put the electrons in their bins based on location
%     x = floor(state(i,1)/1e-9);
%     y = floor(state(i,2)/1e-9);
%     
%     if(x == 0)
%         x = 1;
%     end
%     if(y == 0)
%         y = 1;
%     end
%     
%     % The velocities are then added to the total count for each bin
%     x_temp_sum(x,y) = x_temp_sum(x,y) + state(i,3)^2;
%     y_temp_sum(x,y) = y_temp_sum(x,y) + state(i,4)^2;
%     num_temp(x,y) = num_temp(x,y) + 1;
% end

% Using the sum of the velocities, the temperature of each bin can be
% calculated
% temp = (x_temp_sum + y_temp_sum).*mn./k./2./num_temp;
% temp(isnan(temp)) = 0;
% temp = temp';
% 
% figure(6)
% x = linspace(1, W/1e-9, (W/1e-9)); 
% y = linspace(1, H/1e-9, (H/1e-9)); 
% pcolor(x,y,temp)
% shading interp
% colormap(jet)
% title('Temperature Map (Pt. 3)');
% xlabel('x (nm)');
% ylabel('y (nm)');

% for i = 1:W
%     for j = 1:H
%         
%     end
% end
% figure

