clear all;
close all;

d=2;
epsil=1e-6;
h=0.3;
p_lim=200;
Klimit=round(50*sqrt(d)/h);

K=[1:Klimit];

r_cutoff=min(sqrt(d),h*sqrt(log(1/epsil)));

for i=1:length(K)
    r(i)=K(i)^(-1/d);
    p(i)=0;
    err(i)=1;
    
    while err(i)>epsil & p(i)<p_lim
        p(i)=p(i)+1;
        
        a(i)=r(i);        
        bstar(i)=(a(i)+sqrt((a(i)^2)+2*p(i)*h*h))/2;
        
        b(i)=min(bstar(i),r(i)+r_cutoff);
        
        err(i)=((2^p(i))/(factorial(p(i)))*((a(i)/h)^p(i))*((b(i)/h)^p(i))*exp(-((a(i)-b(i))^2)/(h*h)));
    end
    
    terms(i)=nchoosek(p(i)-1+d,d);

    n(i)=min(K(i),ceil((r_cutoff/r(i))^d));
    C(i)=K(i)+log(K(i))+((1+n(i))*terms(i));
 
    disp(sprintf('d=%d r=%f K=%d r=%f n=%f p=%d terms=%d C=%d ',d,r_cutoff,K(i),r(i),n(i),p(i),terms(i),C(i)));
end

figure;
subplot(2,2,1);
plot(K,r,'r-');
xlabel('K');
ylabel('Maximum cluster radius');

subplot(2,2,2);
plot(K,p,'r-');
xlabel('K');
ylabel('Truncation number');

subplot(2,2,3);
semilogy(K,terms,'r-');hold on
semilogy(K,K,'b-');
xlabel('K');
ylabel('# of terms');

subplot(2,2,4);
semilogy(K,C,'r-');
hold on;
%plot(K,C1,'r--');
xlabel('K');
ylabel('The constant term in the complexity');
[val,index]=min(C);hold on;
v=axis;
X=[K(index) K(index)];
Y=[v(3) v(4)];
line(X,Y);
plot(K(index),v(3),'x');
title(sprintf('d=%d h=%2.2f K=%d p=%d r=%f',d,h,K(index),p(index),r(index)));
    
figure;
fig=plot(K,C,'k-');set(fig,'linewidth',2);
hold on;
fig=xlabel('K');set(fig,'fontsize',20);
fig=ylabel('c');set(fig,'fontsize',20);
fig=gca;
set(fig,'fontsize',15);
box on; grid on;    
filename=sprintf('K_strategy.eps');
print('-depsc','-r300',filename);
    