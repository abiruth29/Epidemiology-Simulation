clear
clc
clearvars
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
susceptible_Humans=zeros(1,no_of_days);
infected_Humans=zeros(1,no_of_days);
Recovered_Humans=zeros(1,no_of_days);
Saa_population = zeros(1,no_of_days);
SAa_population = zeros(1,no_of_days);
SAA_population = zeros(1,no_of_days);
Iaa_population = zeros(1,no_of_days);
IAa_population = zeros(1,no_of_days);
IAA_population = zeros(1,no_of_days);
Naa_population = zeros(1,no_of_days);
NAa_population = zeros(1,no_of_days);
NAA_population = zeros(1,no_of_days);
total_Susceptible_vectors=zeros(1,no_of_days);
total_Infected_vectors=zeros(1,no_of_days);
total_Human_population=zeros(1,no_of_days);
total_Vector_Population=zeros(1,no_of_days);
frequency_aa=zeros(1,no_of_days);
frequency_Aa=zeros(1,no_of_days);
frequency_AA=zeros(1,no_of_days);
frequency_p=zeros(1,no_of_days);
frequency_q=zeros(1,no_of_days);
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
    


    susceptible_Humans(1,i)=Sh0;
    infected_Humans(1,i)=Ih0;
    Recovered_Humans(1,i)=Rh0;
    Saa_population(1,i)=Saa0;
    SAa_population(1,i)=SAa0;
    SAA_population(1,i)=SAA0;
    Iaa_population(1,i)=Iaa0;
    IAa_population(1,i)=IAa0;
    IAA_population(1,i)=IAA0;
    Naa_population(1,i)=Naa0;
    NAa_population(1,i)=NAa0;
    NAA_population(1,i)=NAA0;
    total_Susceptible_vectors(1,i)=Sv0;
    total_Infected_vectors(1,i)=Iv0;
    total_Human_population(1,i)=Nh0;
    total_Vector_Population(1,i)=Nv0;
    frequency_aa(1,i)=Faa0;
    frequency_Aa(1,i)=FAa0;
    frequency_AA(1,i)=FAA0;
    frequency_p(1,i)=p;
    frequency_q(1,i)=q;
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

% simulation graphs
%simulation 1
figure(1)
clf;
plot(Time,susceptible_Humans,LineWidth=2,DisplayName='susceptible');
hold on
plot(Time,infected_Humans,LineWidth=2,DisplayName='infected');
hold on
plot(Time,Recovered_Humans,LineWidth=2,DisplayName='recovered');
hold on
plot(Time,total_Human_population,LineWidth=2,DisplayName='total');
legend('show')
title('humans vs time',FontSize=14);
grid on
hold off

figure(2)
clf;
plot(Time,total_Susceptible_vectors,LineWidth=2,DisplayName='susceptible');
hold on
plot(Time,total_Infected_vectors,LineWidth=2,DisplayName='infected');
hold on
plot(Time,total_Vector_Population,LineWidth=2,DisplayName='total');
legend('show')
title('vectors vs time',FontSize=14);
grid on

figure(3)
clf;
plot(Time,Saa_population,LineWidth=2,DisplayName='susceptible aa');
hold on
plot(Time,Iaa_population,LineWidth=2,DisplayName='infected aa');
hold on
plot(Time,Naa_population,LineWidth=2,DisplayName='total aa');
legend('show')
title('aa allele vectors',FontSize=14);
grid on

figure(4)
clf;
plot(Time,SAa_population,LineWidth=2,DisplayName='susceptible Aa');
hold on
plot(Time,IAa_population,LineWidth=2,DisplayName='infected Aa');
hold on
plot(Time,NAa_population,LineWidth=2,DisplayName='total Aa');
legend('show')
title('Aa allele vectors',FontSize=14);
grid on

figure(5)
clf;
plot(Time,SAA_population,LineWidth=2,DisplayName='susceptible AA');
hold on
plot(Time,IAA_population,LineWidth=2,DisplayName='infected AA');
hold on
plot(Time,NAA_population,LineWidth=2,DisplayName='total AA');
legend('show')
title('AA allele vectors',FontSize=14);
grid on

figure(6)
clf;
plot(Time,Saa_population,LineWidth=2,DisplayName='susceptible aa');
hold on
plot(Time,SAa_population,LineWidth=2,DisplayName='susceptible Aa');
hold on
plot(Time,SAA_population,LineWidth=2,DisplayName='susceptible AA');
legend('show')
title('susceptible vectors vs time',FontSize=14);
grid on

figure(7)
clf;
plot(Time,Iaa_population,LineWidth=2,DisplayName='infected aa');
hold on
plot(Time,IAa_population,LineWidth=2,DisplayName='infected Aa');
hold on
plot(Time,IAA_population,LineWidth=2,DisplayName='infected AA');
legend('show')
title('infected vectors vs time',FontSize=14);
grid on


figure(8)
clf;
plot(Time,frequency_aa,LineWidth=2,DisplayName='frequency aa');
hold on
plot(Time,frequency_Aa,LineWidth=2,DisplayName='frequency Aa');
hold on
plot(Time,frequency_AA,LineWidth=2,DisplayName='frequency AA');
legend('show')
title('Frequency of genotype',FontSize=14);
grid on

figure(9)
clf;
plot(Time,frequency_p,LineWidth=2,DisplayName='frequency p');
hold on
plot(Time,frequency_q,LineWidth=2,DisplayName='frequency q');
legend('show')
title('dominant recessive allele frequency',FontSize=14);
grid on


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

