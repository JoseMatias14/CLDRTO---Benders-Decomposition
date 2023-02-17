
function [Kalman] = FKalman(ny,A,C,it)

V=.5;
W=.5;
sM=size(A);
PP=eye(sM(1));
VV=eye(ny)*V; 
WW=eye(sM(1))*W;
for j=1:it;
    PP = A*PP*A'-A*PP*C'*inv(VV+C*PP*C')*C*PP*A'+ WW;
end
Kalman = A*PP*C'*inv(VV+C*PP*C');



