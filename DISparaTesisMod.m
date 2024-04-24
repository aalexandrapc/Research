function [ddISsB,dIS,dISx,dISy,ddIS,dISMX,dISMY,r,c,dxIS,dyIS,nfd,ncd]=DISparaTesisMod                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%     Densidad construída a partir de la imagen satelital
load('Datacn1.mat');                                                       % this files contains the initial dataset of the bathymetry H                  
[H,Hx,Hy,domx,domy]=batimetria_valencia_Saa(Datacn1);
[NDVIT,DXIS,DYIS,nrN,ncN,BwT,Bww1,Bww2]= ICmodParaTesisDoc ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NDVIT(~(NDVIT>0.01 & NDVIT<1))=0;                                          %NDVIT(~(NDVIT>0.2 & NDVIT<0.5))=0;
[r,c,v]=find(NDVIT(2:end-1,2:end-1));
armis=DXIS*DYIS;                                                           %area de cada rectangulo en la malla de la imagen satelital
TNDVIT=NDVIT+(5*ones(size(NDVIT)));                                        %traslacion NDVIT en dos  unidades hacia arriba
TNDVIT(~BwT)=0;                                                            %se hace cero todolo que está fuera del lago
[m,n]=size(TNDVIT);
N=zeros(1,m*n);
for j=1:m*n
    N(j)=armis*TNDVIT(j);
end
cnd=sum(N) ;                                                               % constante de normalización de densidad
dIS=(1/cnd)*(TNDVIT);                                                      % densidad generada a partir de la imagen satelital                                                           
ddIS=[r,c,v];                                                              %(r,c) son la fila y columna con indice de clorofila entre 0.01 y 1 y v guarda el valor delíndicede clorofila
dIS=double(dIS);                                                           % fue necesaria esta conversión para poderla graficar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Gradiente de la densidad construida a partir de la imagen satelital
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dISx=zeros(m-2,n-2);                                                       % matriz que será llenada con los valores de la derivada de dIS con respecto a x
dISy=zeros(m-2,n-2);                                                       % matriz que será llenada con los valores de la derivada de dIS con respectoa y
for i=2:(m-1)
    for j=2:(n-1)
        dISx(i-1,j-1)=(dIS(i,j+1)-dIS(i,j-1))/(2*DXIS);                    %calculo de la derivada con respecto a x                   
    end
end

for j=2:(n-1)
    for i=2:(m-1)
        dISy(i-1,j-1)=(dIS(i+1,j)-dIS(i-1,j))/(2*DYIS);                    % cálculo de la derivada con respecto a y                      
    end
end
dIS=dIS(2:end-1,2:end-1); 
NDVITT=NDVIT(2:end-1,2:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Puntos con indice de clorofila mayor a 0 no ubicados en borde del lago
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
for i=1:length(ddIS(:,1))
    celda=ddIS(i,1:2);
    aux=NDVITT(celda(1)-1:celda(1)+1,celda(2)-1:celda(2)+1);
    if (length(find(aux))>=7)
        h=h+1;
        ddISsB(h,1:3)=ddIS(i,1:3);%ddISsB(h,1:2)=ddIS(i,1:2);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Densidades marginales construidas a partir de la imagen satelital
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dISMX=sum(DYIS*dIS);                                                       %Densidad marginal de X a partir de dIS
dISMY=sum(DXIS*dIS,2);                                                     %Densidad marginal de Y a partir de dIS
[nfd,ncd]=size(dIS);                                                       % número de filas y columnas de dIS respectivamente 
dxIS=(DXIS:DXIS:DXIS*ncd);                                                 %length(dIS(1,:)));  coordenadas x correspondientes a la malla de laimagen satelital
dyIS=(DYIS:DYIS:DYIS*nfd);                                                 %length(dIS(:,1)) coordenadas y correspondientes a la malla de laimagen satelital                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gráficas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                         


%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[XI,YI]=meshgrid(ddISsB(:,1),ddISsB(:,2));
ZI=griddata(ddISsB(:,1),ddISsB(:,2),ddISsB(:,3),XI,YI);
figure(5500) 
%contour(dxIS,dyIS,dIS,'k');
%hold on
mesh(ZI);
shading interp
figure(5501) 
%contour(dxIS,dyIS,dIS,'k');
%hold on
surf(ddISsB(:,1),ddISsB(:,2),ZI);
shading interp
title('Puntos sobre el lago donde se registro clorofila sin incluir el borde del lago')
hold off
%meshgrid(dxIS,dyIS,dIS)
% figure(3)
% plot(dxIS,dISMX)  
% title('densidad marginal con respecto a X')
% figure(4)
% plot(dyIS,dISMY)
% title('densidad marginal con respecto a Y')
% figure (13), title('derivada con respecto a x')
% surf(dxIS,dyIS,dISx)
% shading interp
% %meshgrid(dxIS,dyIS,dISx)
% figure (14) , title('derivada con respecto a y')
% surf(dxIS,dyIS,dISy)
% shading interp
% %meshgrid(dxIS,dyIS,dISy)

 end