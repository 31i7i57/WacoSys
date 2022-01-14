%% ******************************************************************
%  _________________________________________________________________
% |             |                                                   |
% |   Name      |  Waco-Texas Resilience Analysis                   |
% |_____________|___________________________________________________|
% |             |  Alex and Jon EELE 490 R                          |
% |             |                                                   |
% |             |                                                   |
% | Description |                                                   |
% |             |                                                   |
% |             |                                                   |
% |_____________|___________________________________________________|
% |             |                                                   |
% |   Author    |  © Maryam Bahramipanah+ ....                      |
% |_____________|___________________________________________________|
% |             |                                                   |
% |  Reference  |                                                   |
% |_____________|___________________________________________________|
% |             |                                                   |
% |Last edition |  30 Nov. 2021                                     |
% |_____________|___________________________________________________|

%% ******************************************************************
%
clear all
close all

%% SYSTEM DATA

% Bus data:
        % Type: 1 - Slack Bus; 2 - PQ Bus; 3 - PV Bus
        % Unknown P/Q are put zero. Voltage guess: 1<0 
        % There are four Voltage levels: 138 kV; 69 kV; 13.2 kV; and 12.5 kV
        % Total spot load active power is 564 MW
        % Loads' reactive powers are selected arbitrarily due to lack of data (1/3 of active power) 
        % There are critical and non-critical loads
        % Two power generation units are considered: 1- PV (5 MW) and 2- Sandy 1080 MW
        % Sandy substation is selected as Slack bus
        % Power base is 1000 MVA
        % Voltage base is 138 kV at the slack 
        % Qmin and Qmax are chosen arbitrarily due to lack of data




%        Bus   Bus    Voltage   Angle     ----Generator-------   --------Load-----
%        No    type   Mag.      Degree     PG(pu)       QG(pu)    PL(pu)   QL(pu)      Qmin       Qmax      
busdata=[ 1     1      1         0         0            0          0        0          -2          2; % Sandy 138 kV
          2     2      1         0         0            0          0        0           0          0; % West 138 kV
          3     2      1         0         0            0          0        0           0          0; % Sanger 138 kV
          4     2      1         0         0            0          0        0           0          0; % South 138 kV
          5     2      1         0         0            0          0        0           0          0; % MM 138 kV
          6     3      1.05      0         0.005        0          0        0           -0.002     0.002; % Robinson 138 kV       
          7     2      1         0         0            0          0        0           0          0; % Hewitt 138 kV
          8     2      1         0         0            0          0        0           0          0; % Woodway 138 kV
          9     2      1         0         0            0          0        0           0          0; % Northwest 138 kV
          10    2      1         0         0            0          0        0           0          0; % North 138 kV
          11    2      1         0         0            0          0        0           0          0; % Colonial 138 kV
          12    2      1         0         0            0          0        0           0          0; % East 138 kV
          13    2      1         0         0            0          0        0           0          0; % Lasalle 138 kV
          14    2      1         0         0            0          0        0           0          0; % Bellmead 138 kV
          15    2      1         0         0            0          0        0           0          0; % North Crest 138 kV
          16    2      1         0         0            0          0        0           0          0; % Elm mot 138 kV
          17    2      1         0         0            0          0        0           0          0; % Baylor 138 kV
          18    2      1         0         0            0          0        0           0          0; % Katy 69 kV
          19    2      1         0         0            0          0        0           0          0; % Elm mot 69 kV
          20    2      1         0         0            0          0.007    0.0023      0          0; % Baylor 13.2 kV Critical and Non-Critical
          21    2      1         0         0            0          0.005    0.0017      0          0; % North 13.2 kV Critical
          22    2      1         0         0            0          0.005    0.0016      0          0; % Colonial 13.2 kV Critical
          23    2      1         0         0            0          0.008    0.0026      0          0; % MM 13.2 kV Critical
          24    2      1         0         0            0          0.001    0.0036      0          0; % North Crest 13.2 kV Critical
          25    2      1         0         0            0          0.008    0.0027      0          0; % North Crest 12.5 kV non-Critical
          26    2      1         0         0            0          0.048    0.016       0          0; % Bellmead 12.5 kV Critical and non-Critical
          27    2      1         0         0            0          0.097    0.032       0          0; % East 12.5 kV Critical and non-Critical
          28    2      1         0         0            0          0.0065   0.0021      0          0; % Katy 12.5 kV Critical and non-Critical
          29    2      1         0         0            0          0.014    0.0046      0          0; % Colonial 12.5 kV non-Critical
          30    2      1         0         0            0          0.0168   0.0056      0          0; % North 12.5 kV non-Critical
          31    2      1         0         0            0          0.0257   0.0085      0          0; % Northwest 12.5 kV Critical and non-Critical
          32    2      1         0         0            0          0.011    0.0036      0          0; % Sanger 12.5 kV Critical and non-Critical
          33    2      1         0         0            0          0.034    0.0115      0          0; % Woodway 12.5 kV non-Critical
          34    2      1         0         0            0          0.034    0.0115      0          0; % MM 12.5 kV non-Critical
          35    2      1         0         0            0          0.034    0.0115      0          0; % Hewitt 12.5 kV non-Critical
          36    2      1         0         0            0          0.055    0.018       0          0; % South 12.5 kV Critical and non-Critical
          37    2      1         0         0            0          0.0493   0.0164      0          0; % West 12.5 kV Critical and non-Critical
          38    2      1         0         0            0          0.006    0.002       0          0]; % Lasalle 12.5 kV Critical and non-Critical


% Line/Transformer Data
        % There is no real-data for lines -- Just the lenghth is known
        % Line parameters are given per length of line
        % The parameters should be modified later



%            Bus    bus     R            X                G            B               AMPACITY      length 
%            nl     nr      p.u.         p.u.             p.u.         p.u.            p.u.           mile  
%                
linedata=[   1      12      0.010014      0.023594        0.0          0.0000014       12             9.03; % 138 kV
             1      2       0.010014      0.023594        0.0          0.0000014       12             8.44; % 138 kV
             1      13      0.010014      0.023594        0.0          0.0000014       12             7.06; % 138 kV
             2      3       0.010014      0.023594        0.0          0.0000014       12             1.28; % 138 kV
             2      4       0.010014      0.023594        0.0          0.0000014       12             3.95; % 138 kV
             2      5       0.010014      0.023594        0.0          0.0000014       12             2.97; % 138 kV
             2      6       0.010014      0.023594        0.0          0.0000014       12             3.55; % 138 kV
             2      12      0.010014      0.023594        0.0          0.0000014       12             11.32; % 138 kV
             3      8       0.010014      0.023594        0.0          0.0000014       12             1.83; % 138 kV
             3      9       0.010014      0.023594        0.0          0.0000014       12             2.37; % 138 kV
             4      5       0.010014      0.023594        0.0          0.0000014       12             5.75; % 138 kV
             4      5       0.010014      0.023594        0.0          0.0000014       12             5.75; % 138 kV
             4      11      0.010014      0.023594        0.0          0.0000014       12             2.36; % 138 kV
             5      6       0.010014      0.023594        0.0          0.0000014       12             4.15; % 138 kV
             5      7       0.010014      0.023594        0.0          0.0000014       12             1.59; % 138 kV
             6      7       0.010014      0.023594        0.0          0.0000014       12             3.11; % 138 kV
             5      7       0.010014      0.023594        0.0          0.0000014       12             1.59; % 138 kV
             9      10      0.010014      0.023594        0.0          0.0000014       12             2.95; % 138 kV
             10     11      0.010014      0.023594        0.0          0.0000014       12             2.44; % 138 kV
             10     12      0.010014      0.023594        0.0          0.0000014       12             4.53; % 138 kV
             12     17      0.010014      0.023594        0.0          0.0000014       12             1.36; % 138 kV
             12     13      0.010014      0.023594        0.0          0.0000014       12             3.12; % 138 kV
             12     14      0.010014      0.023594        0.0          0.0000014       12             3.89; % 138 kV
             14     15      0.010014      0.023594        0.0          0.0000014       12             3.39; % 138 kV
             14     16      0.010014      0.023594        0.0          0.0000014       12             4.08; % 138 kV
             18     19      0.010014      0.023594        0.0          0.0000014       8              4.08; % 69 kV
             2      37      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             3      32      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             4      36      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             5      34      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             5      23      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/13.2 kV Transformer
             7      35      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             8      33      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             9      31      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             10     30      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             10     21      0.014546      0.014747        0.0          0.00000059      6              1.00; % 138/13.2 kV Transformer
             11     29      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             11     22      0.014546      0.014747        0.0          0.00000059      6              1.00; % 138/13.2 kV Transformer
             12     18      0.030911      0.031337        0.0          0.00000111      6              1.00; % 138/69 kV Transformer
             12     27      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             13     38      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             14     26      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             15     25      0.023638      0.023964        0.0          0.00000096      6              1.00; % 138/12.5 kV Transformer
             15     24      0.014546      0.014747        0.0          0.00000059      6              1.00; % 138/13.2 kV Transformer
             16     19      0.030911      0.031337        0.0          0.00000111      6              1.00; % 138/69 kV Transformer
             17     20      0.014546      0.014747        0.0          0.00000059      6              1.00; % 138/13.2 kV Transformer
             18     28      0.014546      0.014747        0.0          0.00000059      6              1.00]; % 138/13.2 kV Transformer


% Number of Buses:
N=length(busdata(:,1));

bus = busdata(:,1);            % Bus Number
type = busdata(:,2);           % Type of Bus 1-Slack, 2-PV, 3-PQ

Sbase=1000;
total_load_P=sum(busdata(:,7))*Sbase;
total_load_Q=sum(busdata(:,8))*Sbase;


%% Hourly load profile --- Creation --- hourly

load('RCIcoeff.mat'); % residential, commercial and Industrial load coefficients

% Enter the share of each residential, commercial and industrial loads for each bus

          %R      C      I     bus#
RCI_Share=[0.00   0.00   0.00;  %1
           0.00   0.00   0.00;  %2
           0.00   0.00   0.00;  %3
           0.00   0.00   0.00;  %4
           0.00   0.00   0.00;  %5
           0.00   0.00   0.00;  %6
           0.00   0.00   0.00;  %7
           0.00   0.00   0.00;  %8
           0.00   0.00   0.00;  %9
           0.00   0.00   0.00;  %10
           0.00   0.00   0.00;  %11
           0.00   0.00   0.00;  %12
           0.00   0.00   0.00;  %13
           0.00   0.00   0.00;  %14
           0.00   0.00   0.00;  %15
           0.00   0.00   0.00;  %16
           0.00   0.00   0.00;  %17
           0.00   0.00   0.00;  %18
           0.00   0.00   0.00;  %19
           0.00   1.00   0.00;  %20
           0.00   1.00   0.00;  %21
           0.00   1.00   0.00;  %22
           0.00   1.00   0.00;  %23
           0.00   1.00   0.00;  %24
           0.44   0.56   0.00;  %25 => 0.44 R  & 0.56 C
           0.05   0.60   0.35;  %26 => 0.05 R  & 0.60 C  & 0.35 I
           0.60   0.40   0.00;  %27 => 0.60 R  & 0.40 C
           0.54   0.46   0.00;  %28 => 0.54 R  & 0.46 C
           1.00   0.00   0.00;  %29 => 1.00 R
           1.00   0.00   0.00;  %30 => 1.00 R
           1.00   0.00   0.00;  %31 => 1.00 R
           0.32   0.68   0.00;  %32 => 0.32 R  & 0.68 C
           0.00   1.00   0.00;  %33
           0.008  0.06   0.932; %34 => 0.008 R & 0.06 R & 0.932 I
           0.00   1.00   0.00;  %35
           0.62   0.62   0.38;  %36 => 0.62 R & 0.38 C
           0.02   0.06   0.92;  %37 => 0.02 R & 0.06 C & 0.92 I
           0.00   1.00   0.00]; %38


% Creating Load Profile - Hourly

for i=1:N
    P_hourly(i,:)= (busdata(i,7)*RCI_Share(i,1)).*Rcoeff + (busdata(i,7)*RCI_Share(i,2)).*Ccoeff + (busdata(i,7)*RCI_Share(i,3)).*Icoeff ;
    Q_hourly(i,:)= (busdata(i,8)*0).*Rcoeff + (busdata(i,8)*0.1).*Ccoeff + (busdata(i,8)*0.9).*Icoeff ;  % not sure if we can assume this
end

%% Yearly load profile
P_yearly={};
Q_yearly={};

load('yearcoeff.mat'); % load yearly coefficient for loads

for i=1:365
    P_yearly{i,1}=P_hourly.*yearcoeff(i,1);  % each cell corresponds to a day in a year (38 bus for 24 hours)
    Q_yearly{i,1}=P_hourly.*yearcoeff(i,1); % each cell corresponds to a day in a year (38 bus for 24 hours)
end

%% Hourly PV profile --- Creation --- hourly

load('PVcoeff_hourly.mat');  % 365 x 24 

PV_yearly=busdata(6,5).*PVcoeff_hourly;  % for bus 6 (the only bus that has a PV)


%% COMPUTE YBUS
Nbr=length(linedata(:,1)); % Number of Branches

j=sqrt(-1); 

Ns = linedata(:,1);
Nr = linedata(:,2);
R = linedata(:,3); 
% R = linedata(:,3).*linedata(:,8);
X = linedata(:,4);
% X = linedata(:,4).*linedata(:,8);
Bc = j*linedata(:,6).*linedata(:,8);

Z = R + j*X; % Total
y= ones(Nbr,1)./Z;      %branch admittance  Y=1/Z

N = max(max(Ns),max(Nr));      % No. of buses
ybus=zeros(N,N); % % Initialise YBus

% Formation of the off-diagonal elements
for k=1:Nbr
    ybus(Ns(k),Nr(k))=ybus(Ns(k),Nr(k))-y(k);
    ybus(Nr(k),Ns(k))=ybus(Ns(k),Nr(k));
end

% formation of the diagonal elements
for  n=1:N
    for k=1:Nbr
        if Ns(k)== n
            ybus(n,n) = ybus(n,n) + y(k) + (Bc(k)/2);
        elseif Nr(k)== n
            ybus(n,n) = ybus(n,n) + y(k) + (Bc(k)/2);
        end
    end
end

G = real(ybus);                % Conductance matrix
B = imag(ybus);                % Susceptance matrix


%% FIND PV /Slack & PQ BUSES
             
PV_Buses = find(type == 3);   % PV and Slack Buses

PQ_Buses = find(type == 2);               % PQ Buses
Npv = length(PV_Buses);                   % No. of PV buses
Npq = length(PQ_Buses);                   % No. of PQ buses

%% DEFINE VARIABLES
Voltage={};

%% MAIN LOOP: LOAD FLOW
for day=1:365 % load flow runs for 365 days
    
    PV_day=PV_yearly(day,:); % PV profile for "day"   1 x 24
    Pload_day=P_yearly{day,1};  % P load for "day" for all buses  38 x 24
    Qload_day=Q_yearly{day,1};  % Q load for "day" for all buses  38 x 24

    clear V_hourly;
    
    for hour=1:24  % for each hour
        
        maxiter=100;
        tolerance = 1; %tolerance
        iteration = 1;
        epsilon=1e-5;
        
        V = busdata(:,3);              % Specified Voltage - Initial Voltage magnitude (guess)
        del = busdata(:,4);            % Voltage Angle - Initial Voltage angle (guess)
        
        Pg = busdata(:,5);             % PG - There is only one PV generation
        Qg = busdata(:,6);             % QG 
        Pg(6,:)=PV_day(1,hour);        % load PV profile for "day & hour"


        Pl = Pload_day(:,hour);             % PLoad for "day & hour"
        QL = Qload_day(:,hour);             % PLoad for "day & hour"
               
        
        Qmin = busdata(:,9);           % Minimum Reactive Power Limit
        Qmax = busdata(:,10);          % Maximum Reactive Power Limit


        Pactual = Pg - Pl;                % P Specified..     38 x 1
        Qactual = Qg - QL;                % Q Specified..     38 x 1

        % NEWTON RAPHSON ALGORITHM (LOAD FLOW)
                % it needs Jacobian function to run

        while (tolerance > epsilon) && (iteration < maxiter)
            fprintf(' iteration = ');
            fprintf('  %8.3f', iteration);
            fprintf('\n');
            P=zeros(N,1);
            Q=zeros(N,1);
            % Calculate P and Q
            for k = 1:N
                for n = 1:N
                    P(k) = P(k) + (V(k)*abs(ybus(k,n))*V(n)*cos(del(k)-del(n)-angle(ybus(k,n))));
                    Q(k) = Q(k) + (V(k)*abs(ybus(k,n))*V(n)*sin(del(k)-del(n)-angle(ybus(k,n))));
                    % OR:
                    % P(k) = P(k) + V(k)* V(n)*(G(k,n)*cos(del(k)-del(n)) + B(k,n)*sin(del(k)-del(n)));
                    % Q(k) = Q(k) + V(k)* V(n)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
                end
            end


            % Checking Qlimit violations
            for n = 1:N
                changeit=0;
                if type(n) == 3   % PV buses
                    QG = Q(n)+QL(n);
                    if QG < Qmin(n)
                        QG = Qmin(n);
                        changeit=1;
                    elseif QG > Qmax(n)
                        QG = Qmax(n);
                        changeit=1;
                    end
                    if changeit==1
                        Qactual(n) = QG - QL(n);
                        type(n) = 2;  %change PV bus to PQ bus
                        PV_Buses = find(type == 3);   % PV and Slack Buses
                        PQ_Buses = find(type == 2);   % PQ Buses
                        Npv = length(PV_Buses);       % No. of PV buses
                        Npq = length(PQ_Buses);       % No. of PQ buses
                    end
                end
            end

            % Compute deltaP and deltaQ:  changes from actual value
            dPa = Pactual-P;
            dQa = Qactual-Q;


            % ----- Remove row corresponding to PV buses and slack from dQa
            % ----- Remove row corresponding to slack from dPa
            dQ = dQa(2:N);
            dQ(PV_Buses) = [];
            dP = dPa(2:N);
            dPQ = [dP; dQ];       % Mismatch Vector


            % Compute Jacobian Matrix
            Jac=Jacobian(V,del,ybus,Npq,PQ_Buses);


            % Compute DeltaV and Deltadel
            DeltaVd = inv(Jac)*dPQ;      % Correction Vector
            ddel = DeltaVd(1:N-1);     % Change in Voltage Angle
            dV = DeltaVd(N:end);       % Change in Voltage Magnitude


            % Updating Voltage Angle and Magnitude
            del(2:N) = ddel + del(2:N);    % Voltage Angle

            h = 1;
            for k = 2:N
                if type(k) == 2  % ONLY FOR PQ buses
                    V(k) = dV(h) + V(k);        % Voltage Magnitude
                    h = h+1;
                end
            end


            iteration = iteration + 1;
            tolerance = max(abs(dPQ));

        end
        
       

        % COMPUTE P & Q from SLACK
        Pslack=0;
        Qslack=0;
        for n = 1:N
            Pslack = Pslack + (V(1)*abs(ybus(1,n))*V(n)*cos(del(1)-del(n)-angle(ybus(1,n))));
            Qslack = Qslack + (V(1)*abs(ybus(1,n))*V(n)*sin(del(1)-del(n)-angle(ybus(1,n))));
        end

        Pslack_actual=Sbase*Pslack;
        Qslack_actual=Sbase*Qslack;

        % Save hourly
        V_hourly(:,hour)=V;
        Pslack_hourly(day,hour)=Pslack;  %pu
        Qslack_hourly(day,hour)=Qslack;  %pu
        Pslack_actual_hourly(day,hour)=Pslack_actual;  %pu
        Qslack_actual_hourly(day,hour)=Qslack_actual;  %pu

    end
   
    Voltage{day,1}=V_hourly;
end

%% Here you need to take a look at Voltage at each hour for all buses and Check if they are between 0.95 to 1.05
%% Pslack_actual_hourly gives you the power that Waco gets from Sandy power plant
%% Until here, we are in normal operation. 
%% Feel free to modify the parameters and load/PV data


%% The next step is to creat one or two different scenarios. For instance double(or more) the load and make the voltages less than 0.95 pu
%% In this case, you need to install DG in your system or other management methods can be considered (your call)





%% ****************************************
% The rest of the code here is for a spot 




%% CALLING BUS DATA
% bus = busdata(:,1);            % Bus Number
% type = busdata(:,2);           % Type of Bus 1-Slack, 2-PV, 3-PQ
% V = busdata(:,3);              % Specified Voltage
% del = busdata(:,4);            % Voltage Angle
% Pg = busdata(:,5);             % PG
% Qg = busdata(:,6);             % QG
% Pl = busdata(:,7);             % PL
% QL = busdata(:,8);             % QL
% Qmin = busdata(:,9);           % Minimum Reactive Power Limit
% Qmax = busdata(:,10);          % Maximum Reactive Power Limit


% Pactual = Pg - Pl;                % P Specified..
% Qactual = Qg - QL;                % Q Specified..


%% NEWTON RAPHSON ALGORITHM (LOAD FLOW)
    % it needs Jacobian function to run
% while (tolerance > epsilon) && (iteration < maxiter)
%     fprintf(' iteration = ');
%     fprintf('  %8.3f', iteration);
%     fprintf('\n');
%     P=zeros(N,1);
%     Q=zeros(N,1);
%     % Calculate P and Q
%     for k = 1:N
%         for n = 1:N
%             P(k) = P(k) + (V(k)*abs(ybus(k,n))*V(n)*cos(del(k)-del(n)-angle(ybus(k,n))));
%             Q(k) = Q(k) + (V(k)*abs(ybus(k,n))*V(n)*sin(del(k)-del(n)-angle(ybus(k,n))));
%             % OR:
%             % P(k) = P(k) + V(k)* V(n)*(G(k,n)*cos(del(k)-del(n)) + B(k,n)*sin(del(k)-del(n)));
%             % Q(k) = Q(k) + V(k)* V(n)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
%         end
%     end
%     
%    
%     % Checking Qlimit violations
%     for n = 1:N
%         changeit=0;
%         if type(n) == 3   % PV buses
%             QG = Q(n)+QL(n);
%             if QG < Qmin(n)
%                 QG = Qmin(n);
%                 changeit=1;
%             elseif QG > Qmax(n)
%                 QG = Qmax(n);
%                 changeit=1;
%             end
%             if changeit==1
%                 Qactual(n) = QG - QL(n);
%                 type(n) = 2;  %change PV bus to PQ bus
%                 PV_Buses = find(type == 3);   % PV and Slack Buses
%                 PQ_Buses = find(type == 2);   % PQ Buses
%                 Npv = length(PV_Buses);       % No. of PV buses
%                 Npq = length(PQ_Buses);       % No. of PQ buses
%             end
%         end
%     end
%    
%     % Compute deltaP and deltaQ:  changes from actual value
%     dPa = Pactual-P;
%     dQa = Qactual-Q;
%    
%     
%     % ----- Remove row corresponding to PV buses and slack from dQa
%     % ----- Remove row corresponding to slack from dPa
%     dQ = dQa(2:N);
%     dQ(PV_Buses) = [];
%     dP = dPa(2:N);
%     dPQ = [dP; dQ];       % Mismatch Vector
%     
%         
%     % Compute Jacobian Matrix
%     Jac=Jacobian(V,del,ybus,Npq,PQ_Buses);
%     
%         
%     % Compute DeltaV and Deltadel
%     DeltaVd = inv(Jac)*dPQ;      % Correction Vector
%     ddel = DeltaVd(1:N-1);     % Change in Voltage Angle
%     dV = DeltaVd(N:end);       % Change in Voltage Magnitude
%     
%     
%     % Updating Voltage Angle and Magnitude
%     del(2:N) = ddel + del(2:N);    % Voltage Angle
%     
%     h = 1;
%     for k = 2:N
%         if type(k) == 2  % ONLY FOR PQ buses
%             V(k) = dV(h) + V(k);        % Voltage Magnitude
%             h = h+1;
%         end
%     end
%     
%         
%     iteration = iteration + 1;
%     tolerance = max(abs(dPQ));
%     
% end


%% COMPUTE P & Q from SLACK
% 
% Pslack=0;
% Qslack=0;
% for n = 1:N
%     Pslack = Pslack + (V(1)*abs(ybus(1,n))*V(n)*cos(del(1)-del(n)-angle(ybus(1,n))));
%     Qslack = Qslack + (V(1)*abs(ybus(1,n))*V(n)*sin(del(1)-del(n)-angle(ybus(1,n))));
% end
% 
% Pslack_actual=Sbase*Pslack;
% Qslack_actual=Sbase*Qslack;



