clc;
clear;
h = 1;  
t = 0:h:100; 
k = 1/10; %how long it takes people to recover --> 1/10 = 10 days
R_naught = 3.5;
b = k*R_naught; %each infected conact infects 2 people a day;
population = 4.9e6;
initial_infections = 127; %https://www.hpsc.ie/a-z/respiratory/coronavirus/novelcoronavirus/casesinireland/epidemiologyofcovid-19inireland/march2020/COVID-19%20Epidemiology%20report%20for%20NPHET%2022.03.2020_v1%20web%20version.pdf 
i_initial = initial_infections/population;
r_rate = b/k; %rate of reproduction

s = zeros(1,length(t));
s(1) = 1;
i = zeros(1,length(t)); %rate of infection
i(1) = i_initial;
r = zeros(1,length(t));
r(1) = 0;

%%start day = 27th march
for n = 2: length(t)
    fprintf("Day %i\n",(n-1));
    if(i(n-1) > 0.015 && i(n-1) <= 0.017)
        fprintf("Restriction measures implemented.\n");
        b = 6.5/20;
        s(n) = s(n-1) - b.*(s(n-1)).*(i(n-1))*h;
        i(n) = i(n-1)+(b.*(s(n-1).*i(n-1)) - k.*(i(n-1))).*h;
        r(n) = r(n-1)+k.*(i(n-1).*h);
    elseif(i(n-1) > 0.017 && i(n-1) <= 0.0185)
        fprintf("Schools close \n");
        b = 6.1/20;
        s(n) = s(n-1) - b.*(s(n-1)).*(i(n-1))*h;
        i(n) = i(n-1)+(b.*(s(n-1).*i(n-1)) - k.*(i(n-1))).*h;
        r(n) = r(n-1)+k.*(i(n-1).*h);
    elseif(i(n-1) > 0.0185 && i(n-1) <= 0.0195)
        fprintf("Level 3 lockdown with facemasks enforced\n");
        b = (7/20*0.67);
        s(n) = s(n-1) - b.*(s(n-1)).*(i(n-1))*h;
        i(n) = i(n-1)+(b.*(s(n-1).*i(n-1)) - k.*(i(n-1))).*h;
        r(n) = r(n-1)+k.*(i(n-1).*h);
    elseif(i(n-1) > 0.0195)
        fprintf("full lockdown \n");
        b = 3.5/20;
        s(n) = s(n-1) - b.*(s(n-1)).*(i(n-1))*h;
        i(n) = i(n-1)+(b.*(s(n-1).*i(n-1)) - k.*(i(n-1))).*h;
        r(n) = r(n-1)+k.*(i(n-1).*h);
    else
        fprintf("Locdown is over.\n");
        b = 1/2;
        s(n) = s(n-1) - b.*(s(n-1)).*(i(n-1))*h;
        i(n) = i(n-1)+(b.*(s(n-1).*i(n-1)) - k.*(i(n-1))).*h;
        r(n) = r(n-1)+k.*(i(n-1).*h);
    end
end
fprintf("*****************************************************");
figure(1); hold on
p1 = plot(t,s); L1 = 's(t)';
p2 = plot(t,i); L2 = 'i(t)';
p3 = plot(t,r); L3 = 'r(t)';
legend([p1;p2;p3],L1,L2,L3);
xlabel('Time (days)');
ylabel('Population %');
title('Covid-19 S,I & R Graph For Ireland ');
dateaxis('x',1, '27-Mar-2020');