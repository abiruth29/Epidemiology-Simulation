clear
clc
clearvars

%simulation 1
%constants 
Nh0=9.5*power(10,3); %initial total human population
Sh0=9*power(10,3); %initial susceptible human population
Ih0=2*power(10,2); %initial infected human population
Rh0=3*power(10,2); %initial recovered human population

c=0.457*power(10,-4); % birth rate of human
Muh=0.457*power(10,-4); %death rate of human
Bh=0.4; %vector to human transmission probability
Y=0.121; % human recovery rate

theta=6.353; % vector oviposition rate
Bv=0.4; %human to vector transmission probability

Sv0=9.5*power(10,3); %initial susceptible vector population
Iv0=9.5*power(10,3); %initial infected vector population
Nv0=1.9*power(10,4); %initial total vector population

Naa0=9.5*power(10,3); %initial aa genotype population of vectors
NAa0=4.75*power(10,3); %initial Aa genotype population of vectors
NAA0=4.75*power(10,3); %initial AA genotype population of vectors

Saa0=4.75*power(10,3); %initial aa genotype susceptible population of vectors
SAa0=2.375*power(10,3); %initial Aa genotype susceptible population of vectors
SAA0=2.375*power(10,3); %initial AA genotype susceptible population of vectors

Iaa0=4.75*power(10,3); %initial aa genotype infected population of vectors
IAa0=2.375*power(10,3);%initial Aa genotype infected population of vectors
IAA0=2.375*power(10,3);%initial AA genotype infected population of vectors
K=2*power(10,5); %carrying capacity of the vectors

Muaa=0.25;%death rate of aa genotype vectors
MuAa=0.25;%death rate of Aa genotype vectors
MuAA=0.25;%death rate of AA genotype vectors

q=(2*Naa0 + NAa0)/2*Nv0; % recessive allele frequency
p=(2*NAA0 + NAa0)/2*Nv0; % dominant allele frequency
Faa0=Naa0/Nv0; %frequency of aa genotype
FAa0=NAa0/Nv0; %frequency of Aa genotype
FAA0=NAA0/Nv0; %frequency of AA genotype
no_of_days=1000; %no of days for the simulation

%arrays to store the values
susceptible_Humans1=zeros(1,no_of_days);
infected_Humans1=zeros(1,no_of_days);
Time=zeros(1,no_of_days);
%Model processing
dt=0.01;       % delta t = 1 day
% y(t)=y(t-1) + dt*(function value at t-1)
for i=1:1:no_of_days 
    Sh=Sh0;
    Ih=Ih0;
    Rh=Rh0;
    Nh=Nh0;
    Saa=Saa0;
    SAa=SAa0;
    SAA=SAA0;
    Iaa=Iaa0;
    IAa=IAa0;
    IAA=IAA0;
    Naa=Naa0;
    NAa=NAa0;
    NAA=NAA0;
    Sv=Sv0;
    Iv=Iv0;
    Nv=Nv0;
    q=(2*Naa0 + NAa0)/(2*Nv0);
    p=(2*NAA0 + NAa0)/(2*Nv0);
    Faa=Faa0;
    FAa=FAa0;
    FAA=FAA0;
    


    susceptible_Humans1(1,i)=Sh0;
    infected_Humans1(1,i)=Ih0;
    Time(1,i)=i;


    Sh0=Sh + dt*fx(c,Nh,Bh,Iv,Sh,Muh);
    Ih0=Ih + dt*gx(Bh,Iv,Nh,Sh,Muh,Ih,Y);
    Rh0=Rh + dt*hx(Y,Ih,Muh,Rh);
    Saa0=Saa + dt*px(q,theta,Nv,K,Bv,Ih,Nh,Saa,Muaa);
    Iaa0=Iaa + dt*qx(Bv,Ih,Nh,Saa,Muaa,Iaa);
    SAa0=SAa + dt*rx(p,q,theta,Nv,K,Bv,Ih,Nh,SAa,MuAa);
    IAa0=IAa + dt*sx(Bv,Ih,Nh,SAa,MuAa,IAa);
    SAA0=SAA + dt*tx(p,theta,Nv,K,Bv,Ih,Nh,SAA,MuAA);
    IAA0=IAA + dt*ux(Bv,Ih,Nh,SAA,MuAA,IAA);
    Sv0=Saa0+SAa0+SAA0;
    Iv0=Iaa0+IAa0+IAA0;
    Naa0=Saa0+Iaa0;
    NAa0=SAa0+IAa0;
    NAA0=SAA0+IAA0;
    Nh0=Sh0+Ih0+Rh0;
    Nv0=Sv0+Iv0;
    Faa0=Naa0/Nv0;
    FAa0=NAa0/Nv0;
    FAA0=NAA0/Nv0;
    
end

%simulation 2
%constants 
Nh0=9.5*power(10,3); %initial total human population
Sh0=9*power(10,3); %initial susceptible human population
Ih0=2*power(10,2); %initial infected human population
Rh0=3*power(10,2); %initial recovered human population

c=0.457*power(10,-4); % birth rate of human
Muh=0.457*power(10,-4); %death rate of human
Bh=0.4; %vector to human transmission probability
Y=0.121; % human recovery rate

theta=6.353; % vector oviposition rate
Bv=0.4; %human to vector transmission probability

Sv0=9.5*power(10,3); %initial susceptible vector population
Iv0=9.5*power(10,3); %initial infected vector population
Nv0=1.9*power(10,4); %initial total vector population

Naa0=9.5*power(10,3); %initial aa genotype population of vectors
NAa0=4.75*power(10,3); %initial Aa genotype population of vectors
NAA0=4.75*power(10,3); %initial AA genotype population of vectors

Saa0=4.75*power(10,3); %initial aa genotype susceptible population of vectors
SAa0=2.375*power(10,3); %initial Aa genotype susceptible population of vectors
SAA0=2.375*power(10,3); %initial AA genotype susceptible population of vectors

Iaa0=4.75*power(10,3); %initial aa genotype infected population of vectors
IAa0=2.375*power(10,3);%initial Aa genotype infected population of vectors
IAA0=2.375*power(10,3);%initial AA genotype infected population of vectors
K=2*power(10,5); %carrying capacity of the vectors

Muaa=0.01;%death rate of aa genotype vectors
MuAa=0.25;%death rate of Aa genotype vectors
MuAA=0.25;%death rate of AA genotype vectors

q=(2*Naa0 + NAa0)/2*Nv0; % recessive allele frequency
p=(2*NAA0 + NAa0)/2*Nv0; % dominant allele frequency
Faa0=Naa0/Nv0; %frequency of aa genotype
FAa0=NAa0/Nv0; %frequency of Aa genotype
FAA0=NAA0/Nv0; %frequency of AA genotype
no_of_days=1000; %no of days for the simulation

%arrays to store the values
susceptible_Humans2=zeros(1,no_of_days);
infected_Humans2=zeros(1,no_of_days);
%Model processing
dt=0.01;       % delta t = 1 day
% y(t)=y(t-1) + dt*(function value at t-1)
for i=1:1:no_of_days 
    Sh=Sh0;
    Ih=Ih0;
    Rh=Rh0;
    Nh=Nh0;
    Saa=Saa0;
    SAa=SAa0;
    SAA=SAA0;
    Iaa=Iaa0;
    IAa=IAa0;
    IAA=IAA0;
    Naa=Naa0;
    NAa=NAa0;
    NAA=NAA0;
    Sv=Sv0;
    Iv=Iv0;
    Nv=Nv0;
    q=(2*Naa0 + NAa0)/(2*Nv0);
    p=(2*NAA0 + NAa0)/(2*Nv0);
    Faa=Faa0;
    FAa=FAa0;
    FAA=FAA0;
    


    susceptible_Humans2(1,i)=Sh0;
    infected_Humans2(1,i)=Ih0;


    Sh0=Sh + dt*fx(c,Nh,Bh,Iv,Sh,Muh);
    Ih0=Ih + dt*gx(Bh,Iv,Nh,Sh,Muh,Ih,Y);
    Rh0=Rh + dt*hx(Y,Ih,Muh,Rh);
    Saa0=Saa + dt*px(q,theta,Nv,K,Bv,Ih,Nh,Saa,Muaa);
    Iaa0=Iaa + dt*qx(Bv,Ih,Nh,Saa,Muaa,Iaa);
    SAa0=SAa + dt*rx(p,q,theta,Nv,K,Bv,Ih,Nh,SAa,MuAa);
    IAa0=IAa + dt*sx(Bv,Ih,Nh,SAa,MuAa,IAa);
    SAA0=SAA + dt*tx(p,theta,Nv,K,Bv,Ih,Nh,SAA,MuAA);
    IAA0=IAA + dt*ux(Bv,Ih,Nh,SAA,MuAA,IAA);
    Sv0=Saa0+SAa0+SAA0;
    Iv0=Iaa0+IAa0+IAA0;
    Naa0=Saa0+Iaa0;
    NAa0=SAa0+IAa0;
    NAA0=SAA0+IAA0;
    Nh0=Sh0+Ih0+Rh0;
    Nv0=Sv0+Iv0;
    Faa0=Naa0/Nv0;
    FAa0=NAa0/Nv0;
    FAA0=NAA0/Nv0;
    
end

%simulation 3
%constants 
Nh0=9.5*power(10,3); %initial total human population
Sh0=9*power(10,3); %initial susceptible human population
Ih0=2*power(10,2); %initial infected human population
Rh0=3*power(10,2); %initial recovered human population

c=0.457*power(10,-4); % birth rate of human
Muh=0.457*power(10,-4); %death rate of human
Bh=0.4; %vector to human transmission probability
Y=0.121; % human recovery rate

theta=6.353; % vector oviposition rate
Bv=0.4; %human to vector transmission probability

Sv0=9.5*power(10,3); %initial susceptible vector population
Iv0=9.5*power(10,3); %initial infected vector population
Nv0=1.9*power(10,4); %initial total vector population

Naa0=9.5*power(10,3); %initial aa genotype population of vectors
NAa0=4.75*power(10,3); %initial Aa genotype population of vectors
NAA0=4.75*power(10,3); %initial AA genotype population of vectors

Saa0=4.75*power(10,3); %initial aa genotype susceptible population of vectors
SAa0=2.375*power(10,3); %initial Aa genotype susceptible population of vectors
SAA0=2.375*power(10,3); %initial AA genotype susceptible population of vectors

Iaa0=4.75*power(10,3); %initial aa genotype infected population of vectors
IAa0=2.375*power(10,3);%initial Aa genotype infected population of vectors
IAA0=2.375*power(10,3);%initial AA genotype infected population of vectors
K=2*power(10,5); %carrying capacity of the vectors

Muaa=0.25;%death rate of aa genotype vectors
MuAa=0.01;%death rate of Aa genotype vectors
MuAA=0.01;%death rate of AA genotype vectors

q=(2*Naa0 + NAa0)/2*Nv0; % recessive allele frequency
p=(2*NAA0 + NAa0)/2*Nv0; % dominant allele frequency
Faa0=Naa0/Nv0; %frequency of aa genotype
FAa0=NAa0/Nv0; %frequency of Aa genotype
FAA0=NAA0/Nv0; %frequency of AA genotype
no_of_days=1000; %no of days for the simulation
%arrays to store the values
susceptible_Humans3=zeros(1,no_of_days);
infected_Humans3=zeros(1,no_of_days);
%Model processing
dt=0.01;       % delta t = 1 day
% y(t)=y(t-1) + dt*(function value at t-1)
for i=1:1:no_of_days 
    Sh=Sh0;
    Ih=Ih0;
    Rh=Rh0;
    Nh=Nh0;
    Saa=Saa0;
    SAa=SAa0;
    SAA=SAA0;
    Iaa=Iaa0;
    IAa=IAa0;
    IAA=IAA0;
    Naa=Naa0;
    NAa=NAa0;
    NAA=NAA0;
    Sv=Sv0;
    Iv=Iv0;
    Nv=Nv0;
    q=(2*Naa0 + NAa0)/(2*Nv0);
    p=(2*NAA0 + NAa0)/(2*Nv0);
    Faa=Faa0;
    FAa=FAa0;
    FAA=FAA0;
    


    susceptible_Humans3(1,i)=Sh0;
    infected_Humans3(1,i)=Ih0;


    Sh0=Sh + dt*fx(c,Nh,Bh,Iv,Sh,Muh);
    Ih0=Ih + dt*gx(Bh,Iv,Nh,Sh,Muh,Ih,Y);
    Rh0=Rh + dt*hx(Y,Ih,Muh,Rh);
    Saa0=Saa + dt*px(q,theta,Nv,K,Bv,Ih,Nh,Saa,Muaa);
    Iaa0=Iaa + dt*qx(Bv,Ih,Nh,Saa,Muaa,Iaa);
    SAa0=SAa + dt*rx(p,q,theta,Nv,K,Bv,Ih,Nh,SAa,MuAa);
    IAa0=IAa + dt*sx(Bv,Ih,Nh,SAa,MuAa,IAa);
    SAA0=SAA + dt*tx(p,theta,Nv,K,Bv,Ih,Nh,SAA,MuAA);
    IAA0=IAA + dt*ux(Bv,Ih,Nh,SAA,MuAA,IAA);
    Sv0=Saa0+SAa0+SAA0;
    Iv0=Iaa0+IAa0+IAA0;
    Naa0=Saa0+Iaa0;
    NAa0=SAa0+IAa0;
    NAA0=SAA0+IAA0;
    Nh0=Sh0+Ih0+Rh0;
    Nv0=Sv0+Iv0;
    Faa0=Naa0/Nv0;
    FAa0=NAa0/Nv0;
    FAA0=NAA0/Nv0;
    
end

% simulation graphs
%simulation 1
figure(1)
clf;
plot(Time,susceptible_Humans1,LineWidth=2,DisplayName='simulation 1');
hold on
plot(Time,susceptible_Humans2,LineWidth=2,DisplayName='simulation 2');
hold on
plot(Time,susceptible_Humans3,LineWidth=2,DisplayName='simulation 3');
legend('show')
title('comparison between susceptible humans in 1,2,3',FontSize=14);
grid on
hold off

figure(2)
clf;
plot(Time,infected_Humans1,LineWidth=2,DisplayName='simulation 1');
hold on
plot(Time,infected_Humans2,LineWidth=2,DisplayName='simulation 2');
hold on
plot(Time,infected_Humans3,LineWidth=2,DisplayName='simulation 3');
legend('show')
title('comparison between infected humans in 1,2,3',FontSize=14);
grid on
hold off


%functions for the model equations
function y1=fx(c,Nh,Bh,Iv,Sh,Muh)
    y1=c*Nh - (Bh*Iv/Nh)*Sh - Muh*Sh;
end

function y2=gx(Bh,Iv,Nh,Sh,Muh,Ih,Y)
    y2=(Bh*Iv/Nh)*Sh - Muh*Ih - Y*Ih;
end

function y3=hx(Y,Ih,Muh,Rh)
    y3=(Y*Ih)-(Muh*Rh);
end

function y4=px(q,theta,Nv,K,Bv,Ih,Nh,Saa,Muaa)
    y4=q*q*theta*Nv*(1 - Nv/K)-(Bv*Ih/Nh)*Saa - Muaa*Saa;
end

function y5=qx(Bv,Ih,Nh,Saa,Muaa,Iaa)
    y5=(Bv*Ih/Nh)*Saa - Muaa*Iaa;
end
function y6=rx(p,q,theta,Nv,K,Bv,Ih,Nh,SAa,MuAa)
    y6=2*p*q*theta*Nv*(1 - Nv/K)-(Bv*Ih/Nh)*SAa - MuAa*SAa;
end
function y7=sx(Bv,Ih,Nh,SAa,MuAa,IAa)
    y7=(Bv*Ih/Nh)*SAa - MuAa*IAa;
end
function y8=tx(p,theta,Nv,K,Bv,Ih,Nh,SAA,MuAA)
    y8=p*p*theta*Nv*(1 - Nv/K)-(Bv*Ih/Nh)*SAA - MuAA*SAA;
end
function y9=ux(Bv,Ih,Nh,SAA,MuAA,IAA)
    y9=(Bv*Ih/Nh)*SAA - MuAA*IAA;
end
