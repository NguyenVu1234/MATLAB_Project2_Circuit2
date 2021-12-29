clear;clc;close all;
% x = 10;
% y = 5;
% u = x > y
% v = x < y
% w = x==y
% ww = x>=2*y
% 
% if x>y
%     y=x
% else
%     x=y
% end
% 
% T = 1.0e-03;
% Vm = 1;
% t = linspace(0,T,1001);  
% timeIncrement = t(2)-t(1)
% % 
% % 
% for ii = 1:length(t)
%      v(ii) = (Vm/T)*t(ii);
% end
% plot(t,v)
% xlabel('t')
% ylabel('v(t)')
% title('v(t) vs t')
% %  
% plot(t,v)
T=10^(-3);
Vm=1;
t=linspace(0,T,1001);
%provides increamental value
tincreament=t(2)-t(1)
% for R=1:length(t)
%     v(R)=(Vm/T)*t(R);
% end
% figure
% plot(t,v);
% xlabel('t');
% ylabel('v(t)');
%     
% 
T = 1.0e-03;
Vm = 1;
t = linspace(-T,2*T,3001);  

for ii = 1:length(t)
    if t(ii) >= 0 && t(ii) <= T
        v(ii) = (Vm/T)*t(ii);
    else 
        v(ii) = 0;
    end
    
end

v = [v(find(v)) v(find(v)) v(find(v)) 0]
plot(t,v)
xlabel('t')
ylabel('v(t)')
title('v(t) vs t')

% amplitude = whatever;
% period = whatever;
% y = mod(x,period);
% ix = y > period/2;
% y(ix) = period - y(ix);
% y = (amplitude * 2 / period) * y;

% T=.001;
% Vm=1;
% t1=linspace(-T,0,3001);
% N=1;
% 
% for R=1:length(t1)
% if t1>=0
% n(R)=mod(t1(R),T);
% v1(R)=Vm-(Vm/T)*n(R);
% else
% n(R)=mod(abs(t1(R)),T);
% v1(R)=(Vm/T)*n(R);
% end
% end
% t2=linspace(0,2*T,3001);
% N=1;
% 
% for R=1:length(t2)
% if t2>=0
% n(R)=mod(t2(R),T);
% v2(R)=Vm-(Vm/T)*n(R);
% else
% n(R)=mod(abs(t2(R)),T);
% v2(R)=(Vm/T)*n(R);
% end
% end
% v=[v1,v2]
% t=[t1,t2]
% plot(t,v)
T=.001;
 Vm=1;
t =linspace(-T,2*T,3001);
% 
% for ii = 1:length(t)
%     if t(ii) >=0 && t(ii) <=T
%         v(ii) =  Vm-(Vm/T)*mod(t(ii),T);
%     else
%         v(ii) = 0;
%     end
% end
% plot(t,v)
 T=.001;
 Vm=1;
t =linspace(-T,2*T,3001);

for ii = 1:length(t)
    if t(ii) >=-T && t(ii) <= 2*T 
        v(ii) = (Vm/T)*mod(t(ii),T);
    else 
        v(ii) = 0;
    end
    
end

plot(t,v)
xlabel('t')
ylabel('v(t)')
title('v(t) vs t')


