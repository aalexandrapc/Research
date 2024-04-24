function [H,Hx,Hy,domx,domy]=batimetria_valencia_Saa(Datacn1);

%%% se debe introdicir el tamaño de los cuadrados %%%

%%%% Dx=1000 para el calculo de la constante

                 %%%% Calculo de los puntos medios %%%%

                 
% ESTOS SON LAS DIMENSIONES DE LA MALLA QUE SE USARA PARA MODELAR EL LAGO.                 
Nx=38;
Ny=30;

% LOS VALORES DE n Y m QUE SE INCORPORAN CORRESPONDEN AL TAMAÑO DE LA CELDA
% QUE SE CONSIDERARA, ES DECIR, DX Y DY.
n=1/3;
m=1/3;
x=0:n:Nx;
y=0:m:Ny;

for i=1:length(x)-1
     xm(i)=(x(i)+x(i+1))/2;
end
for i=1:length(y)-1
     ym(i)=(y(i)+y(i+1))/2;
end
 
                 %%%% Calcula la batimetría %%%%
                    
% PARA LA BATIMETRIA, LOS VALORES DE Z DEBEN SER NEGATIVOS. DE ESTA FORMA
% SE REQUIEREN EN EL CALCULO DE LAS VELOCIDADES.

x=Datacn1(:,1);
y=Datacn1(:,2);
z=-Datacn1(:,3); 

[XI,YI]=meshgrid(xm,ym);

ZI=griddata(x,y,z,XI,YI);


for i=1:1:length(ZI(:,1))
    xin=ZI(i,:);
    xin(isnan(xin))=4.5;
    ZI(i,:)=xin;        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:length(ZI(:,1))
    for j=1:1:length(ZI(1,:))
        if ZI(i,j) == 0;
            ZI(i,j)=4.5;
        end
    end
end
%ZI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=10:110                                                               %j=11:1:109
    for i=16:83                                                            %i=17:1:82
        Hm(i-15,j-9)=ZI(i,j);                                              %Hm(i-16,j-10)=ZI(i,j);
    end
end
%Hm
Hx=zeros(length(Hm(:,1))-2,length(Hm(1,:))-2);
Hy=zeros(length(Hm(:,1))-2,length(Hm(1,:))-2);
for i=2:67
    for j=2:100
        Hx(i-1,j-1)=(Hm(i,j+1)-Hm(i,j-1))/(2*n);                           %Hx(i-1,j-1)=(Hm(i,j+1)-Hm(i,j-1))/(2*n);
    end
end

for j=2:100
    for i=2:67
        Hy(i-1,j-1)=(Hm(i+1,j)-Hm(i-1,j))/(2*m);                           %Hy(i-1,j-1)=(Hm(i+1,j)-Hm(i-1,j))/(2*m);
    end
end
H=Hm(2:(length(Hm(:,1))-1),2:(length(Hm(1,:))-1));
domx=(n:n:length(H(1,:))*n);
domy=(m:m:length(H(:,1))*m);
figure(25)
 surf(domx,domy,H);
%  figure(26)
%  [d1,d2]=meshgrid(domx,domy)
%  mesh(d1,d2,H);
%   figure(26)
%   surf(domy,domx,H');
 