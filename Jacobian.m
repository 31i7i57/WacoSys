%% ******************************************************************
%  _________________________________________________________________
% |             |                                                   |
% |   Name      |  Jacobian matrix calculation                      |
% |_____________|___________________________________________________|
% |             |  This is a general code for Jacobian matrix       |
% |             |  calculation for Power flow calculations          |
% |             |  Using NEWTON-RAPHSON method                      |
% | Description |  course: EELE 454- Power system Analysis          |
% |             |                                                   |
% |             |                                                   |
% |_____________|___________________________________________________|
% |             |                                                   |
% |   Author    |  Maryam Bahramipanah                              |
% |_____________|___________________________________________________|
% |             |                                                   |
% |  Reference  |  Power system analysis and design by Glover       |
% |_____________|___________________________________________________|
% |             |                                                   |
% |Last edition |  15 Nov. 2019                                     |
% |_____________|___________________________________________________|

%% ******************************************************************
%
function [Jac]=Jacobian(V,del,ybus,Npq,PQ_Buses)
N=length(ybus);

% COMPUTE J1 - Derivative of Real Power Injections with Angles  
J1 = zeros(N-1,N-1);
for row = 1:(N-1)
    k = row+1;
    for coloumn = 1:(N-1)
        n = coloumn+1;
        if n == k   %J1 -Diagonal n=k
            for n = 1:N
%               J1(row,coloumn) = J1(row,coloumn) + -V(k)* V(n)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
                J1(row,coloumn) = J1(row,coloumn) + (-V(k)*V(n)*abs(ybus(k,n))*sin(del(k,1)-del(n,1)-angle(ybus(k,n))));
            end
%             J1(row,coloumn) = J1(row,coloumn) - V(k)^2*B(k,k);
            J1(row,coloumn) = J1(row,coloumn) + V(k)^2*abs(ybus(k,k))*sin(-angle(ybus(k,k)));
        else   %J1 -off Diagonal n~=k
%           J1(row,coloumn) = V(k)* V(n)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
            J1(row,coloumn) = V(k)* V(n)* abs(ybus(k,n))*sin(del(k)-del(n)-angle(ybus(k,n)));
        end
    end
end

 
% COMPUTE J2 - Derivative of Real Power Injections with V
J2 = zeros(N-1,Npq);
for row = 1:(N-1)
    k = row+1;
    for coloumn = 1:Npq
        n = PQ_Buses(coloumn);
        if n == k %J2 -Diagonal n=k
            for n = 1:N
%                 J2(row,coloumn) = J2(row,coloumn) + V(n)*(G(k,n)*cos(del(k)-del(n)) + B(k,n)*sin(del(k)-del(n)));
                J2(row,coloumn) = J2(row,coloumn) + V(n)*abs(ybus(k,n))*cos(del(k)-del(n)-angle(ybus(k,n)));
            end
%             J2(row,coloumn) = J2(row,coloumn) + V(k)*G(k,k);
            J2(row,coloumn) = J2(row,coloumn) + V(k)*abs(ybus(k,k))*cos(angle(ybus(k,k)));
        else %J2 -off Diagonal n~=k
%             J2(row,coloumn) = V(k)*(G(k,n)*cos(del(k)-del(n)) + B(k,n)*sin(del(k)-del(n)));
            J2(row,coloumn) = V(k,1)*abs(ybus(k,n))*cos(del(k)-del(n)-angle(ybus(k,n)));
        end
    end
end



% COMPUTE J3 -Derivative of Reactive Power Injections with Angles
J3 = zeros(Npq,N-1);
for row = 1:Npq
    k = PQ_Buses(row);
    for coloumn = 1:(N-1)
        n = coloumn+1;
        if n == k
            for n = 1:N
%                 J3(row,coloumn) = J3(row,coloumn) + V(k)* V(n)*(G(k,n)*cos(del(k)-del(n)) + B(k,n)*sin(del(k)-del(n)));
                J3(row,coloumn) = J3(row,coloumn) + (V(k)*V(n)*abs(ybus(k,n))*cos(del(k,1)-del(n,1)-angle(ybus(k,n))));
            end
%             J3(row,coloumn) = J3(row,coloumn) - V(k)^2*G(k,k);
            J3(row,coloumn) = J3(row,coloumn) - V(k)^2*abs(ybus(k,k))*cos(-angle(ybus(k,k)));
        else
%             J3(row,coloumn) = V(k)* V(n)*(-G(k,n)*cos(del(k)-del(n)) - B(k,n)*sin(del(k)-del(n)));
            J3(row,coloumn) = -V(k)* V(n)* abs(ybus(k,n))*cos(del(k)-del(n)-angle(ybus(k,n)));
        end
    end
end


% COMPUTE J4 - Derivative of Reactive Power Injections with V
J4 = zeros(Npq,Npq);
for row = 1:Npq
    k = PQ_Buses(row);
    for coloumn = 1:Npq
        n = PQ_Buses(coloumn);
        if n == k %J4 -Diagonal n=k
            for n = 1:N
%                 J4(row,coloumn) = J4(row,coloumn) + V(n)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
                J4(row,coloumn) = J4(row,coloumn) + V(n)*abs(ybus(k,n))*sin(del(k)-del(n)-angle(ybus(k,n)));
            end
%             J4(row,coloumn) = J4(row,coloumn) - V(k)*B(k,k);
            J4(row,coloumn) = J4(row,coloumn) - V(k)*abs(ybus(k,k))*sin(angle(ybus(k,k)));
        else %J4 -off Diagonal n~=k
%             J4(row,coloumn) = V(k)*(G(k,n)*sin(del(k)-del(n)) - B(k,n)*cos(del(k)-del(n)));
            J4(row,coloumn) = V(k,1)*abs(ybus(k,n))*sin(del(k)-del(n)-angle(ybus(k,n)));
        end
    end
end


% Jacobian Matrix
Jac=[J1 J2;J3 J4];













