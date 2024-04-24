function [DIXLHM1P,DIYLHM1P,bvnd]=DILHM1P(cbH,fbH,cb,fb,H,DXIS,DYIS,dISconBorde, dIS,dISx,dISy,PIXLHM,PIYLHM,L);
%% DILHM1P determina la posicion en la malla batimetria_valencia_Saa de puntos 
%% cuya posición original esta determinada en función de la malla de la  IS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nfd,ncd]=size(dIS);                                                       %número de filas y columnas de dIS respectivamente 
dxIS=(DXIS:DXIS:DXIS*ncd);                                                 %length(dIS(1,:)));  coordenadas x correspondientes a la malla de laimagen satelital
dyIS=(DYIS:DYIS:DYIS*nfd);                                                 %length(dIS(:,1)) coordenadas y correspondientes a la malla de laimagen satelital                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Gráfica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                         
figure (1)
contour(dxIS,dyIS,dISconBorde,[0 0],'black');
hold on
contour(dxIS,dyIS,dIS,[0,0],'black');
for i=1:length(PIXLHM)
    plot(PIXLHM(i),PIYLHM(i),'*')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% La función que viene genera las trayectorias del LHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [XELHM,YELHM, XEm, YEm,contLHM,L,bvnd]=LangevinHM(DXIS,DYIS,dIS,dISx,dISy,PIXLHM,PIYLHM,L); 

%% Grafico de trayectorias
figure(2)
contour(dxIS,dyIS,dIS,[0,0],'black'); 
hold on
for i=1:length(XELHM(:,1))
    plot(XELHM(i,:),YELHM(i,:))
    plot(XELHM(i,1),YELHM(i,1),'.','MarkerEdgeColor','r')
     plot(XELHM(i,end),YELHM(i,end),'.','MarkerEdgeColor','black')
end
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cambio de malla
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cont=contLHM;
DX=1/3;
DY=1/3;
distX=XEm'-DXIS*(repmat(cb(1),[cont,1]));                                        %Distancia entre el extremo izquierdo del lago y un punto con clorofila  
distY=YEm'-DYIS*(repmat(fb(1),[cont,1]));                                        %Distancia entre el extremo i del lanferior del lago y un punto con clorofila  
DIXLHM1P=DX*(repmat(cbH(1),[cont,1]))+distX;
DIYLHM1P=DY*(repmat(fbH(1),[cont,1]))+distY;

%% Gráficas
domx=(DX:DX:length(H(1,:))*DX);                                            
domy=(DY:DY:length(H(:,1))*DY);
figure(3)
 contour(domx,domy,H,[0 0],'black')
  hold on
  contour(domx,domy,H,4.5,'red');
 for i=1:1:length(DIXLHM1P)
   plot(DIXLHM1P(i),DIYLHM1P(i))
 end
hold off
end
