function [Sn] = Roads_Add_HBEFA_Parameters(S)
%--------------------------------------------------------------------------
% Function to add the HBEFA road properties needed to assign emission
% factors to different roads. This conversion is partly based on the road
% type and size, partly on the speed limit and for some road types, we use
% the ADT to determine whether it is a larger or smaller type of road.
%
% 20.09.2018 -Henrik Grythe
% Kjeller NILU 
%--------------------------------------------------------------------------

global HBEFA_roads
fprintf('------------------------------------------------------------------\n')
fprintf('* in Add_Road_HBEFA_Parameters  *\n')
fprintf('---   Reading available classifications \n')

fprintf('Convertion based on convertion sheet:\n%s\n',HBEFA_roads)
[Rl,~,RlT] = xlsread(HBEFA_roads,'Convert_RTMroads_to_HBEFA');

Rname = RlT(2:end,1);

% Accepted Input for VK is: all else is "null" 
rl=[{'E'},{'R'},{'F'},{'K'},{'P'}];

% extract the necessary fields
URB = extractfield(S,'URBAN');
SPD = extractfield(S,'SPEED');
DIS = extractfield(S,'DISTANCE');
ADT = extractfield(S,'SUM_ADT');
RDT = extractfield(S,'VK');

RDI=zeros(size(RDT));
for i=1:length(rl)
    I=find(ismember(RDT,rl(i)));
    RDI(I)=i;
end
% Unclassified roads are here treated as "Municipal Roads" RDI = 4. 
RDI(RDI==0)=4;

% Workaround for the 5 roads that have an erranous speed in the input Roads recieved from Urbanet on 2018/06/27. 
% Theres a bug in the files with  unrealistic and unassigned speeds.
SPD         = round(SPD/10)*10;
SPD(SPD<30) = 30;
SPD(SPD>110)=110;


% counter (used to keep track so that an easy check if all roads have been classified can be done)
t=0; 

% Loop all roads to see what matching class can be made for HBEFA
% 1st Test: URBAN / RURAL
% 2nd Test: SIGNED SPEED LIMIT
% 3rd Test: Type of Road.
% 4th Test: ADT
% 
HBEFA_NUM  = zeros(size(SPD));
HBEFA_EQIV = zeros(size(SPD));
for i=1:length(URB)
    if URB(i) % )Access, )Distr, )Local, )MW-City, )MW-Nat.,)Trunk-City, )Trunk-Nat
        Cn=0;
        if     SPD(i)==30 || SPD(i)==40 %1)Access road
            N=1;
        elseif SPD(i)==50               %)Access ) Distr )Local )Trunk-City
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
            
        elseif SPD(i)==60 % 2)Distr, 3)Local, 4)MW-City, 6)Trunk-City
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
            
        elseif SPD(i)==70   % 2)Distr, 4)MW-City, 6)Trunk-City, 7)Trunk-Nat
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
        elseif SPD(i)==80   % 2)Distr, 4)MW-City, 5)MW-Nat., 6)Trunk-City, 7)Trunk-Nat
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=3; end
            end
        elseif SPD(i)==90   % 4)MW-City, 5)MW-Nat., 6)Trunk-City, 7)Trunk-Nat
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
        elseif SPD(i)==100 ||  SPD(i)==110 % 4)MW-City, 5)MW-Nat., 7)Trunk-Nat
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=3; end
            end
        else
            error(sprintf('Speed not matching speedlimits defined SPD %i ',SPD(i)))
        end
    else % IS RURAL  % 1)Access, 2)Distr, 3)Distr sin., 4)Local, 5)Local sin., 6)Trunk, 7)Semi-MW, 8)MW
        Cn=0;
        if     SPD(i)==30 || SPD(i)==40 %1)Access road
            N=1;
        elseif SPD(i)==50               % 1)Access, 2)Distr, 3)Distr sin., 4)Local, 5)Local sin.
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=3; end
            end
        elseif SPD(i)==60 || SPD(i)==70 % 2)Distr, 3)Distr sin., 4)Local, 5)Local sin. 6)Trunk
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=3; end
            end
        elseif SPD(i)==80               % 2)Distr, 3)Distr sin., 4)Local, 5)Local sin. 6)Trunk 8)MW
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
        elseif SPD(i)==90               % 2)Distr, 3)Distr sin.,  6)Trunk 7)Semi-MW  8)MW
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
        elseif SPD(i)==100              % 2)Distr, 3)Distr sin.,  6)Trunk   8)MW
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=4; end
            end
        elseif SPD(i)==110              % 4)MW-City, 5)MW-Nat., 7)Trunk-Nat
            if RDI(i)==1 || RDI(i)==2 || RDI(i)==3
                N=1;
            elseif RDI(i)==4 || RDI(i)==5
                if ADT(i)>600; N=2; elseif ADT(i)<=600 && ADT(i)>300; N=3; else; N=3; end
            end
        else
            error('Speed not matching speedlimits defined')
        end

    end
    % Test that all roads are present
    ix=find(Rl(:,5)==SPD(i) & Rl(:,2)==URB(i)& Rl(:,7)==N & Rl(:,6)==Cn);
    if length(ix)~=1
        fprintf('SPD: %i URB:%i N:%i antall: %i \n',SPD(i),URB(i),N, length(ix))
        HBEFA_NUM(i)  = Rl(ix(1),1);
        HBEFA_EQIV(i) = Rl(ix(1),3);
    else
        HBEFA_NUM(i)  = Rl(ix,1);
        HBEFA_EQIV(i) = Rl(ix,3);
    end
    if rem(i,100000)==0 ; fprintf('Processed: %i of %i RoadLinks\n',i,length(SPD)); end
end % loop links
fprintf('Adding fields to Table ...\n')
T = struct2table(S);
T.HBEFA_NUM  = HBEFA_NUM';
T.HBEFA_EQIV = HBEFA_EQIV';
fprintf('Structuring Table ...\n')
Sn = table2struct(T);
fprintf('Setfields\n  -- HBEFA_NUM\n  -- HBEFA_EQIV ...on all links \n')
fprintf('* end HEDGE_Add_Road_HBEFA_Parameters  *\n')
fprintf('------------------------------------------------------------------\n')
end

