%
% record of revisions:
%  15/07/2017
%  
%
% D = sqrt( (1.2-xr)^2+(0.5-yr)^2+1.87^2);
% cos(phi) = 1.87/sqrt( (1.2-xr)^2+(0.5-yr)^2+1.87^2)=1.87/D
%  
% m = -log(2)/log(cos(phi/2))
% h =((m+1)*Ar/2*pi*D*D)*cos(phi)^m*cos(psi)= ((m+1)*Ar/2*pi*D*D)*(1.87/D)^m*(f/(1.15*D));
% f = 0.543*(1.2-xr)+0.198*(0.5-yr)+1.87
% mod(vn) = 1.15;
% 1.87^2= 3.497
 %[x,y] = meshgrid(0:0.5:5, 0:0.5:3);

 clear all;clc;
%  xr, yr set to 101 points respectively;
x = 0:0.05:5;
y = 0:0.03:3;
[xr,yr] = meshgrid(x,y);
 
 % PD effective area
Ar = 1.5e-3;

 % calculate the distance between the LED1 and (xr,yr),that is D1
 % Xrt1 = [ a1,b1,1.87]
 % For LED2,LED3,LED4,the same as the LED1
a1 = 1.2-xr;
a2 = 3.8-xr;
a3 = 1.2-xr;
a4 = 3.8-xr;

b1 = 0.5-yr;
b2 = 0.5-yr;
b3 = 2.5-yr;
b4 = 2.5-yr;

c1 = a1.^2;
d1 = b1.^2;
c2 = a2.^2;
d2 = b2.^2;
c3 = a3.^2;
d3 = b3.^2;
c4 = a4.^2;
d4 = b4.^2;

D1 = sqrt(c1+d1+3.497);
D2 = sqrt(c2+d2+3.497);
D3 = sqrt(c3+d3+3.497);
D4 = sqrt(c4+d4+3.497);


 % vdp --vector dot product,that is Vrt*Vn   
 % Vrt = (x-xr,y-yr,z-zr)=(a-xr,b-yr,1.87)
 % Vn = (0.543,0.198,1)
vdp1 = 0.543*a1+0.198*b1+1.87;
vdp2 = 0.543*a2+0.198*b2+1.87;
vdp3 = 0.543*a3+0.198*b3+1.87;
vdp4 = 0.543*a4+0.198*b4+1.87;

% calculate the Lambertian parameters of four LEDs,that is m1,m2,m3,m4
% cos(phi/2) = sqrt( ( 1+cos(phi) )/2=sqrt((1+1.87./D)/2)
 
m1 =  -log(2)./log(sqrt((1+1.87./D1)/2));
m2 =  -log(2)./log(sqrt((1+1.87./D2)/2));
m3 =  -log(2)./log(sqrt((1+1.87./D3)/2));
m4 =  -log(2)./log(sqrt((1+1.87./D4)/2));


% calculate the LOS link gain h for four LEDs,
% cos(psi) =(Vrt*Vn)/(|Vrt|*|Vn|) =vdp./(1.15*D)
% cos(phi) = |zt-zr|/sqrt((xt-xr).^2+(yt-yr).^2+(zt-zr).^2)=1.87./D
h1 = (Ar*(m1+1)./(2*pi*D1.^2)).*(vdp1./(1.15*D1)).*(1.87./D1).^m1;
h2 = (Ar*(m2+1)./(2*pi*D2.^2)).*(vdp2./(1.15*D2)).*(1.87./D2).^m2;
h3 = (Ar*(m3+1)./(2*pi*D3.^2)).*(vdp3./(1.15*D3)).*(1.87./D3).^m3;
h4 = (Ar*(m4+1)./(2*pi*D4.^2)).*(vdp4./(1.15*D4)).*(1.87./D4).^m4;

%%%%%%%%%%%%%%
% maximum of h1,h2,h3,h4 is the link gain and the other three as the
% interference

 
h11 = zeros(101,101);
h12 = zeros(101,101);
h13 = zeros(101,101);
h14 = zeros(101,101);
h21 = zeros(101,101);
h22 = zeros(101,101);
h23 = zeros(101,101);
h24 = zeros(101,101);
h31 = zeros(101,101);
h32 = zeros(101,101);
h33 = zeros(101,101);
h34 = zeros(101,101);
h41 = zeros(101,101);
h42 = zeros(101,101);
h43 = zeros(101,101);
h44 = zeros(101,101);

%to calculate the largest link gain for 4 LEDs, the other three link gain
%as interference
% h11 is the LOS link gain for LED1 , and the other h12,h13,h14 as the
% interference
for ii = 1:101
    for jj = 1:101
        
if h1(ii,jj)>=h2(ii,jj)&&h1(ii,jj)>=h3(ii,jj)&&h1(ii,jj)>=h4(ii,jj)
    h11(ii,jj) = h1(ii,jj);
    h12(ii,jj) = h2(ii,jj);
    h13(ii,jj) = h3(ii,jj);
    h14(ii,jj) = h4(ii,jj);
elseif h2(ii,jj)>=h1(ii,jj)&&h2(ii,jj)>=h3(ii,jj)&&h2(ii,jj)>=h4(ii,jj)
    h21(ii,jj) = h1(ii,jj);
    h22(ii,jj) = h2(ii,jj);
    h23(ii,jj) = h3(ii,jj);
    h24(ii,jj) = h4(ii,jj);
elseif h3(ii,jj)>=h1(ii,jj)&&h3(ii,jj)>=h2(ii,jj)&&h3(ii,jj)>=h4(ii,jj)
    h31(ii,jj) = h1(ii,jj);
    h32(ii,jj) = h2(ii,jj);
    h33(ii,jj) = h3(ii,jj);
    h34(ii,jj) = h4(ii,jj);
elseif h4(ii,jj)>=h1(ii,jj)&&h4(ii,jj)>=h2(ii,jj)&&h4(ii,jj)>=h3(ii,jj)
    h41(ii,jj) = h1(ii,jj);
    h42(ii,jj) = h2(ii,jj);
    h43(ii,jj) = h3(ii,jj);
    h44(ii,jj) = h4(ii,jj);
else
end

    end
end
 
%%% calculate the sum of interference for four LEDs
% h1j = h12+h13+h14 --the total interference for LED1
% similarily, h2j,h3j,h4j
h1j = zeros(101,101);
h2j = zeros(101,101);
h3j = zeros(101,101);
h4j = zeros(101,101);
for ii = 1:101;
    for jj = 1:101;
    if  h11(ii,jj)~=0
        h1j(ii,jj) = h12(ii,jj)+h13(ii,jj)+h14(ii,jj);
    else
    end
    end
end
for ii = 1:101;
    for jj = 1:101;
    if  h22(ii,jj)~=0
        h2j(ii,jj) = h21(ii,jj)+h23(ii,jj)+h24(ii,jj);
    else
    end
    end
end

for ii = 1:101;
    for jj = 1:101;
    if  h33(ii,jj)~=0
        h3j(ii,jj) = h31(ii,jj)+h32(ii,jj)+h34(ii,jj);
    else
    end
    end
end
for ii = 1:101;
    for jj = 1:101;
    if  h44(ii,jj)~=0
        h4j(ii,jj) = h41(ii,jj)+h42(ii,jj)+h43(ii,jj);
    else
    end
    end
end
 

q = 1.6e-19;
B =1e6;
ksi = 0.5;   % decision threshold
Rp = 2.933e-2;
p = 51.3*2;
pj = p*4;
Ibg = 1e-8;
I2 = 0.562;


sigmaI1 =  0.25*(pj^2)*(Rp^2)*h1j+2*q*Ibg*I2*B;
sigmaI2 =  0.25*(pj^2)*(Rp^2)*h2j+2*q*Ibg*I2*B;
sigmaI3 =  0.25*(pj^2)*(Rp^2)*h3j+2*q*Ibg*I2*B;
sigmaI4 =  0.25*(pj^2)*(Rp^2)*h4j+2*q*Ibg*I2*B;
ksai1 =  0.5*pj*Rp*h1j;
ksai2 =  0.5*pj*Rp*h2j;
ksai3 =  0.5*pj*Rp*h3j;
ksai4 =  0.5*pj*Rp*h4j;
kkk = ksai1;
 


%
z10 = zeros(101,101);
for ii = 1:101
    for jj = 1:101
      if sigmaI1(ii,jj)>=1e-6
         z10(ii,jj) =  (ksi - ksai1(ii,jj))./sqrt(sigmaI1(ii,jj));
      else  
        z10(ii,jj) = 0;
      end
    end 
end

z11 =zeros(101,101);
for ii = 1:101
    for jj = 1:101
        if sigmaI1(ii,jj)>=1e-6
            z11(ii,jj) = (p*Rp*h11(ii,jj)+ksai1(ii,jj)-ksi)./sqrt(sigmaI1(ii,jj)+2*q*p*Rp*B*h11(ii,jj));
        else  
            z11(ii,jj) = 0;
        end
    end
end

% test2 = z10;
% test1 = z11;
%        
% mmmm = 0.5*qfunc(z10);
% nnn = 0.5*qfunc(z11);


z20 = zeros(101,101);
for ii = 1:101
    for jj = 1:101
      if sigmaI2(ii,jj)>=1e-6
         z20(ii,jj) =  (ksi - ksai2(ii,jj))./sqrt(sigmaI2(ii,jj));
      else 
        z20(ii,jj) = 0;
      end
    end 
end

z21 =zeros(101,101);
for ii = 1:101
    for jj = 1:101
        if sigmaI2(ii,jj)>=1e-6
            z21(ii,jj) = (p*Rp*h22(ii,jj)+ksai2(ii,jj)-ksi)./sqrt(sigmaI2(ii,jj)+2*q*p*Rp*B*h22(ii,jj));
        else  
            z21(ii,jj) = 0;
        end
    end
end

z30 = zeros(101,101);
for ii = 1:101
    for jj = 1:101
      if sigmaI3(ii,jj)>=1e-6
         z30(ii,jj) =  (ksi - ksai3(ii,jj))./sqrt(sigmaI3(ii,jj));
      else  
        z30(ii,jj) = 0;
      end
    end 
end

z31 =zeros(101,101);
for ii = 1:101
    for jj = 1:101
        if sigmaI3(ii,jj)>=1e-6
            z31(ii,jj) = (p*Rp*h33(ii,jj)+ksai3(ii,jj)-ksi)./sqrt(sigmaI3(ii,jj)+2*q*p*Rp*B*h33(ii,jj));
        else  
            z31(ii,jj) = 0;
        end
    end
end



z40 = zeros(101,101);
for ii = 1:101
    for jj = 1:101
      if sigmaI4(ii,jj)>=1e-6
         z40(ii,jj) =  (ksi - ksai4(ii,jj))./sqrt(sigmaI4(ii,jj));
      else 
        z40(ii,jj) = 0;
      end
    end 
end

z41 =zeros(101,101);
for ii = 1:101
    for jj = 1:101
        if sigmaI4(ii,jj)>=1e-6
            z41(ii,jj) = (p*Rp*h44(ii,jj)+ksai4(ii,jj)-ksi)./sqrt(sigmaI4(ii,jj)+2*q*p*Rp*B*h44(ii,jj));
        else  
            z41(ii,jj) = 0;
        end
    end
end
 
ber11 = 0.5*qfunc(z10)+0.5*qfunc(z11);
ber22 = 0.5*qfunc(z20)+0.5*qfunc(z21);
ber33 = 0.5*qfunc(z30)+0.5*qfunc(z31);
ber44 = 0.5*qfunc(z40)+0.5*qfunc(z41);
ber = ber11+ber22+ber33+ber44;

v = [0.01,0.001,1e-5];
out11 = contour(xr,yr,ber,3);
% hold on
% out22 = contour(xr,yr,ber22,[3]);
% hold on
% out33 = contour(xr,yr,ber33,[3]);
% hold on
% out44 = contour(xr,yr,ber44,[3]);

%[0.01,0.001,1e-5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure(2)  intensity;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
L = 6000;
mm1 = ones(101,101);
 
LL1 = (mm1+1);
 
I1 = (L*LL1.*((1.87./D1).^LL1))./(2*pi*D1.^2);
I2 = (L*LL1.*((1.87./D2).^LL1))./(2*pi*D2.^2);
I3 = (L*LL1.*((1.87./D3).^LL1))./(2*pi*D3.^2);
I4 = (L*LL1.*((1.87./D4).^LL1))./(2*pi*D4.^2);
I = I1+I2+I3+I4;
 figure(2)
surf(xr,yr,I);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 
