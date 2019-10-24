function [I,J,x0,y0] = findMatrixIndex(Ori,xy)
% Ori <- <struct>;
% Ori.LT
% Ori.RB
% Ori.dx
% Ori.dy

X = xy(1);
Y = xy(2);
I = round((Ori.LT(2)-Y)/Ori.dy)+1;
J = round((X-Ori.LT(1))/Ori.dx)+1;

x0 = (J-1)*Ori.dx+Ori.LT(1);
y0 = Ori.LT(2)-(I-1)*Ori.dy;
end