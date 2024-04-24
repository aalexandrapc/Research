function [W1,W2]=MBT(L,NT,DT)

dW = sqrt(DT)*randn(NT,2*L);
dW1 = dW(:,1:L);
dW2 = dW(:,L+1:2*L);
W1 =  [zeros(NT,1),cumsum(dW1,2)];
W2 =  [zeros(NT,1),cumsum(dW2,2)];