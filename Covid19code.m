NumberOfSims = 30;
SusceptibleSim = cell(NumberOfSims,1);
InfectedSim = cell(NumberOfSims,1);
RecoveredSim = cell(NumberOfSims,1);

for sim = 1:NumberOfSims
    
    limits_of_Area = 100;         
    People = 100;         
    Underlying_Conditions = People*0.04;% 4% of people with underlying health conditions would need to be
    %hospitilized - https://www.thelancet.com/coronavirus/archive
    p_underlying_health_cond = 0.2;
    Health_Cond_Test = People*p_underlying_health_cond;
    Health_NonCond_Test = People*(1-p_underlying_health_cond);
    
    People_size = 0.95;            % Size of People
    speed = 10;         % People speed
    t_recovery = 10;    % Mean recovery time (in days)
    days = 30;
    Time_step = 0.02;        % Time increment (in days)
    isolation_rate = 0.1;          % % of people in isolation
    initial_infection_rate = 0.1;       % Percent of carriers initially infected
    transmission_rate_strain1 = 0.50;     % Transmission rte of Covid Strain
    %transmission_rate_strain2 = 0.99;
    mortality_rate_strain1 = 0.03;      % Death rate
    %mortality_rate_strain2 = 0.025;
    mortality_rate_health_conditions = 0.1078; %Underlying condition mortality : based on Ireland Data
    positions = rand(People,2).*(limits_of_Area);                 % Initial Positions of People
    isolated = (rand(People,1)<isolation_rate);              % People in isolation
    Initial_Direction = rand(People,1).*2.*pi;                  % Direction of carriers
    People_Movement_Speed = [speed.*cos(Initial_Direction) speed.*sin(Initial_Direction)];         % Velocity of people
    
    infected = rand(People,1)<initial_infection_rate;
    healthy = ~infected;                    %Healthy people
    healthy_health_conditions = rand(People,1)< p_underlying_health_cond;
    recovered = zeros(People,1);                 % Recovered people
    dead = zeros(People,1);                      % Deaths
    p_death = rand(People,1)<mortality_rate_strain1; % Probability that a person will die
    p_death_Health_conditions = rand(Health_Cond_Test,1)<mortality_rate_health_conditions; %Probability that a person WITH HEALTH CONDITIONS will die
    p_death_Health_conditions2 = zeros(People,1);
    p_death_Health_conditions(numel(p_death_Health_conditions2)) = 0;
    
    
    t_rec = ceil(t_recovery.*ones(People,1)+4.*randn(People,1))./Time_step; % Time-to-recover from COVID - Randomized
    t = linspace(0,days,days./Time_step);      % Time
    collision = zeros(People,People);    % Matrix to keep track of recent collisions
    
    
    % Allocate space for solution vectors
    inf_sum = zeros(days/Time_step,0);          % % of infected carriers
    hea_sum = zeros(days/Time_step,0);          % % of healthy carriers
    rec_sum = zeros(days/Time_step,0);          % % of recovered carriers
    dead_sum = zeros(days/Time_step,0);         % e of dead carriers
    rem_sum = zeros(days/Time_step,0);
    cumulative_sum = zeros(days/Time_step,0);   % % of cumulative disease case
    for i = 1:(days/Time_step)
        % Update carrier position
        pos_new_of_person = positions+People_Movement_Speed.*(~repmat(isolated,1,2)).*Time_step;
        
        % Step through each carrier
        for k = 1:People
            
            % If recovery time is up, carrier is either recovered or dead
            if infected(k)&&t_rec(k)<=0
                
                % If recovery time is up and carrier is dead, well, it's dead.
                % Zero it's velocity
                if p_death(k) || p_death_Health_conditions(k)
                    dead(k) = 1;
                    People_Movement_Speed(k,:) = [0 0];
                    recovered(k)=0;
                    infected(k)=0;
                    healthy(k)=0;
               
                else
                    recovered(k)=1;
                    infected(k)=0;
                    healthy(k)=0;
                    dead(k)=0;
                end
                
                % If carrier is infected and not recovered, decrement recovery time
            elseif (infected(k))
                t_rec(k) = t_rec(k)-1;
            end
            
            % Step through all other carriers, looking for collisions, and if
            % so, transmit disease and recalculate trajectory
            for j = 1:1:People
                if j~=k
                    % Get positions of carriers j and k
                    pos1 = pos_new_of_person(k,:);
                    pos2 = pos_new_of_person(j,:);
                    
                    % If collision between two living specimens, re-calcuate
                    % direction and transmit virus (but don't check the same
                    % two carriers twice)
                    if norm(pos1-pos2)<=(2*People_size) && ~collision(k,j) && ~collision(j,k)
                        phi = atan2((pos2(2)-pos1(2)),(pos2(1)-pos1(1)));
                           
                        if isolated(k)||dead(k)
                            
                            % Get normal direction vector of 'virtual wall'
                            phi_wall = -phi+pi/2;
                            n_wall = [sin(phi_wall) cos(phi_wall)];
                            dot = People_Movement_Speed(j,:)*n_wall';
                        else
                            
                            % Get velocity magnitudes
                            %v1_mag = sqrt(People_Movement_Speed(k,1)^2+People_Movement_Speed(k,2)^2);
                            %v2_mag = sqrt(People_Movement_Speed(j,1)^2+People_Movement_Speed(j,2)^2);
                            
                            % Get directions
                            %th1 = atan2(People_Movement_Speed(k,2),People_Movement_Speed(k,1));
                            %th2 = atan2(People_Movement_Speed(j,2),People_Movement_Speed(j,1));
                            
                        end
                        
                        % If either is infected and not dead...
                        if((infected(j)||infected(k)))&&((~dead(k)||~dead(j)))
                            
                            % If either is recovered, no transmission
                            if recovered(k) || isolated(k) && healthy(k)
                                infected(k)=0;
                            elseif recovered(j) || isolated(j) && healthy(j)
                                infected(j)=0;
                               
                                % Otherwise, transmit virus
                            else
                                    Covid_Data = rand(1)<transmission_rate_strain1;
                                    if Covid_Data
                                        infected(j)=1;
                                        infected(k)=1;
                                        healthy(j)=0;
                                        healthy(k)=0;
                                    end
%                                 else 
%                                     Covid_Data = rand(1) < transmission_rate_strain2;
                                %end
                            end
                        end
                    end
                end
            end
            
            % Look for collisions with outer walls and re-direct
            
            % Left Wall
            if pos_new_of_person(k,1)<=People_size
                if People_Movement_Speed(k,1)<0
                    People_Movement_Speed(k,1)=-People_Movement_Speed(k,1);
                end
                
                % Right wall
            elseif pos_new_of_person(k,1)>=limits_of_Area-People_size
                if People_Movement_Speed(k,1)>0
                    People_Movement_Speed(k,1)=-People_Movement_Speed(k,1);
                end
            end
            
            % Bottom Wall
            if pos_new_of_person(k,2) <=People_size
                if People_Movement_Speed(k,2)<0
                    People_Movement_Speed(k,2)=-People_Movement_Speed(k,2);
                end
                
                % Top Wall
            elseif pos_new_of_person(k,2) >=limits_of_Area-People_size
                if People_Movement_Speed(k,2)>0
                    People_Movement_Speed(k,2)=-People_Movement_Speed(k,2);
                end
            end
            
        end
        
        
        
        
        % Update color vector
        color = [infected healthy recovered].*(1-dead);
        
        % Update solution vectors
        inf_sum(i) = sum(infected)*100/People;
        hea_sum(i) = sum(healthy)*100/People;
        rec_sum(i) = sum(recovered)*100/People;
        dead_sum(i) = sum(dead)*100/People;
        cumulative_sum(i) = 100-hea_sum(i);
        
        
        % Initialize plots on first loop iteration
        if i==1
            % Plot transmission simulation
            figure(1);
            %lineWidth = 2;
            markerSize = 16;
            
            subplot(2,3,[1 2 4 5]);
            h = scatter(pos_new_of_person(:,1),pos_new_of_person(:,2),markerSize,color,'filled','MarkerEdgeColor','k'); hold on;
            
            xlim([0,limits_of_Area]);
            ylim([0,limits_of_Area]);
            axis square;
            grid on;
            box on;
            set(gca,'YTickLabel',[]);
            set(gca,'XTickLabel',[]);
            titlestring = strcat('Percent: Susceptible= ',num2str(hea_sum(i)),', Infected= ',...
                num2str(inf_sum(i)),', Recovered=', num2str(rec_sum(i)),...
                ', Deceased=', num2str(dead_sum(i)));
            title(titlestring);
            
            % Resize markers to match infection radius
            currentunits = get(gca,'Units');
            set(gca, 'Units', 'Points');
            axpos = get(gca,'Position');
            set(gca, 'Units', currentunits);
            markerWidth = 2*People_size/diff(xlim)*axpos(3); % Calculate Marker width in points
            lineWidth = 0.5*markerWidth;
            set(h, 'SizeData', markerWidth^2);
            markerSize = markerWidth^2;
            
            % Plot infection rates vs. time
            subplot(2,3,[3 6]);
            h2 = plot(t(1:i),hea_sum(1:i),'g','LineWidth',2);hold on;
            h3 = plot(t(1:i),inf_sum(1:i),'r','LineWidth',2);hold on;
            h4 = plot(t(1:i),rec_sum(1:i),'b','LineWidth',2);hold on;
            h5 = plot(t(1:i),dead_sum(1:i),'k','LineWidth',2);hold on;
            % h6 = plot(t(1:i),dead_health_conditions_sum(1:i),'m','LineWidth',2);hold off;
            legend('Susceptible','Infected','Recovered','Deceased');%,'Deceased-PoorHealth');
            xlabel('Days');
            ylabel('Percent of Population');
            xlim([0,days]);
            ylim([0,100]);
            set(gcf,'Color','w');
            
            % Update data on subesequent iterations
        else
            subplot(2,3,[1 2 4 5]);
            set(h,'XData',pos_new_of_person(:,1));
            set(h,'YData',pos_new_of_person(:,2));
            set(h,'CData',color);
            
            
            % Update title
            titlestring = strcat('Percent: Susceptible= ',num2str(hea_sum(i)),', Infected= ',...
                num2str(inf_sum(i)),', Recovered=', num2str(rec_sum(i)),...
                ', Deceased=', num2str(dead_sum(i)));%,...
            %', Deceased with Health Conditions =' ,num2str( dead_health_conditions_sum(i)));
            title(titlestring);
            
            subplot(2,3,[3 6]);
            set(h2,'XData',t(1:i)); set(h2,'YData',hea_sum(1:i));
            set(h3,'XData',t(1:i)); set(h3,'YData',inf_sum(1:i));
            set(h4,'XData',t(1:i)); set(h4,'YData',rec_sum(1:i));
            set(h5,'XData',t(1:i)); set(h5,'YData',dead_sum(1:i));
            % set(h6,'XData',t(1:i)); set(h6,'YData',dead_health_conditions_sum(1:i));
            
        end
        drawnow;
        positions = pos_new_of_person;
    end
    % Save results to a cell array
    Covid_Data = {'time',t};
    Covid_Data(2,:) = {'unaffected',hea_sum};
    Covid_Data(3,:) = {'infected',inf_sum};
    Covid_Data(4,:) = {'recovered',rec_sum};
    Covid_Data(5,:) = {'deceased',dead_sum};
    % Covid_Data(6,:) = {'deceased-poorhealth',dead_health_conditions_sum};
    SusceptibleSim{sim,1} = hea_sum;
    InfectedSim{sim,1} = inf_sum;
    RecoveredSim{sim,1} = rec_sum+dead_sum;
    hea_sum = 0;
    inf_sum = 0;
    rec_sum = 0;
    dead_sum = 0;
    close all;
end
SusceptibleArray = zeros(length(t),1);
InfectedArray = zeros(length(t),1);
RemovedArray = zeros(length(t),1);

SusSum = 0;
InfSum = 0;
RemSum = 0;
for size1 = 1 : length(t)
    for size = 1:NumberOfSims
        sus_sum = SusceptibleSim{size,1};
        inf_sum = InfectedSim{size,1};
        rem_sum = RecoveredSim{size,1};
        
        SusSum = SusSum + sus_sum(size1);
        InfSum = InfSum + inf_sum(size1);
        RemSum = RemSum + rem_sum(size1);
    end
    SusceptibleArray(size1,1) = SusSum/NumberOfSims;
    InfectedArray(size1,1) = InfSum/NumberOfSims;
    RemovedArray(size1,1) = RemSum/NumberOfSims;
    
    SusSum = 0;
    InfSum = 0;
    RemSum = 0;
end
plot(t,SusceptibleArray,t,InfectedArray,t,RemovedArray)
xlabel('Time (Days) ');
ylabel('% of Population');
title('SIR Simulator Graph');
legend('Susceptible','Infected','Removed');
