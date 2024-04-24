function [XELHM,YELHM, XEm, YEm,contLHM,L,bvnd]=LangevinHM(DXIS,DYIS,dIS,dISx,dISy,PIXLHM,PIYLHM,L);                   
%% dISy,dISx:            Derivadas con respecto a x e y de la densidad dIS construída a partir de la imagen satelital
%% W1 y W2:              Movimientos  brownianos independientes
%% ddIS:                 Matriz de 3 columas que guarda los pares ordenados con índice de clorofila mayor a cero y sus valores
%% L                     Número de pasos de tiempo asociados a la trayectoria   
%% PIXTDBM y PIYTDBM     Estos vectores guardan las coordenadas x e y de los puntos sobre los conjuntos con 0<dIS<1, donde inician las trayectorias de la aproximación de  Langevin-Hasting-Metropolis
[IMAX,JMAX]=size(dIS);                                                     %% IMAX contiene el número de filas de dIS y JMAX su número de columnas
NT =length(PIXLHM);                                                        %% Número de trayectorias a simular                                                                     
DT = 5/10000;                                                              %% Intervalo de tiempo por cada paso de la trayectoria
[W1,W2]=MBT(L,NT,DT);                                                      %% W1 y W2 son aproximaciones de   movimientos Brawniones independientes unidimensionales de media 0 y varianza DT, cada uno
%% Estimación de Trayectorias 
DX=DXIS;
DY=DYIS;
XE = zeros(NT,L+1);                                                        %% Esta matríz contendrá las coordenadas x de la  'cadena candidata' 
YE = zeros(NT,L+1);                                                        %% Esta matríz contendrá las coordenadas y de la  'cadena candidata' 
XE(:,1)=PIXLHM;                                                            %% Esta matríz contendrá las coordenadas x de la aproximación usando Hasting-Metropolis
YE(:,1)=PIYLHM;                                                            %% Esta matríz contendrá las coordenadas y de la aproximación  usando Hasting-Metropolis 
XELHM=zeros(NT,L);
YELHM=zeros(NT,L);
XELHM(:,1)=PIXLHM;
YELHM(:,1)=PIYLHM; 
dISLHM= zeros(NT,L);                                                       %% Esta matriz guardará los valores de dIS en los puntos correspondientes a la discretizacion de la ecuación  de Langevin ajustada por el HM
drifLHMx= zeros(NT,L);                                                     %% Esta matriz guardará los valores de dISx en los puntos correspondientes a la discretizacion de la ecuación  de Langevin ajustada por el HM
drifLHMy= zeros(NT,L);                                                     %% Esta matriz guardará los valores de dISy en los puntos correspondientes a la discretizacion de la ecuación  de Langevin ajustada por el HM                          
for i=1:1:NT
    for j=1:1:L
        JCELT = length((DX:DX:XE(i,j)));
        ICELT = length((DY:DY:YE(i,j)));
      
        PARTX = (XE(i,j)-JCELT*DX)/DX;
        PARTY = (YE(i,j)-ICELT*DY)/DY;
        
        %% A CONTINUACIÓN SE COLOCA LA FUNCIÓN QUE CALCULA LA INTERPOLACIÓN DE U, V Y H
        [dISxPOS,dISyPOS,dISPOS]=interpolacionesDBM(dISx,dISy,dIS,ICELT,JCELT,PARTX,PARTY); 
        dISLHM(i,j)=dISPOS;                                                                                                                            %% Densidad invariante evaluada en (XE(i,j),YE(i,j)
        drifLHMx(i,j)=(1/2)*(dISxPOS/dISPOS);
        drifLHMy(i,j)=(1/2)*(dISyPOS/dISPOS);                                                                                                           %% driff evaluado en (XE(i,j),YE(i,j)
        if (j>=2)
            dTChainf=normpdf(XE(i,j),XELHM(i,j-1)+(drifLHMx(i,j-1)*DT),(DT)^(1/2))*normpdf(YE(i,j),YELHM(i,j-1)+(drifLHMy(i,j-1)*DT),(DT)^(1/2));       %% densidad de transición de la cadena candidata
            dTChainb=normpdf(XELHM(i,j-1),XE(i,j)+(drifLHMx(i,j)*DT),(DT)^(1/2))*normpdf(YELHM(i,j-1),YE(i,j)+(drifLHMy(i,j)*DT),(DT)^(1/2));                                                    
            FSEL=(dISLHM(i,j)*dTChainb)/(dISLHM(i,j-1)*dTChainf);                                                                                  %% Factor de selección
            PSEL=min(1,FSEL);   
            PSELU=rand(1);                                                                                                                         %% Probabilidad de selección 
            if (PSELU < PSEL)
                XELHM(i,j)=XE(i,j);
                YELHM(i,j)=YE(i,j);
            else        
                XELHM(i,j)=XELHM(i,j-1);
                YELHM(i,j)=YELHM(i,j-1);
                drifLHMx(i,j)= drifLHMx(i,j-1);
                drifLHMy(i,j)= drifLHMy(i,j-1);
            end
        end
        XPOS=XELHM(i,j) +(drifLHMx(i,j)*DT) + (W1(i,j+1)-W1(i,j));                                     
        YPOS=YELHM(i,j) + (drifLHMy(i,j)*DT) + (W2(i,j+1)-W2(i,j));                                
        
        %% A CONTINUACIÓN SE COLOCA LA FUNCIÓN QUE VERIFICA SI UNA PARTICULA ESTA EN AGUA O EN TIERRA
         if ((XPOS >= DX) && (XPOS <= DX*JMAX) && (YPOS >= DY) && (YPOS <= DY*IMAX))  
            TRACER=[XELHM(i,j),YELHM(i,j)];
            [XE,YE]=ReboteSimetricoDBM(i,j,dIS,DX,DY,XPOS,YPOS,TRACER,XE,YE);                            
         else                                                                             
            XE(i,j+1)=XELHM(i,j);
            YE(i,j+1)=YELHM(i,j);
        end   
        
    end
end
bvnd=length(find(isnan([XELHM,YELHM])));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Puntos seleccionados de 7 trayectorias 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contLHM=0;
 for i=1:NT
     for j=1:100:L
        contLHM=contLHM+1;
        XEm(contLHM)=XELHM(i,j);
        YEm(contLHM)=YELHM(i,j);
     end
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A CONTINUACIÓN SE COLOCA LA FUNCIÓN QUE CALCULA LA DENSIDAD POR NUCLEO

%[dK,xd,yd]=densidad_nucleo(H,L,XE,YE,NT,cx,cy);

% A CONTINUACIÓN SE COLOCA LAS GRAFICAS PARA CADA CALCULO
% [dmx,dmy]=meshgrid(dxIS,dyIS);
% [dmx,dmy]=meshgrid(dxIS,dyIS);
% figure(1)
% hold on
% % for i=1:length(r)
%     plot(DX*c',DY*r,'r');
% % end
% figure(2)
% plot(XE(:,100+1),YE(:,100+1),'.','MarkerEdgeColor','r')
% % hold on
% % contour(dmx,dmy,dIS,'k');
% % hold off
%  figure(3)
% hold on
% plot(XE(:,500+1),YE(:,500+1),'.','MarkerEdgeColor','r')
% % hold on
% % contour(dmx,dmy,dIS,'k');
% % hold off
%  figure(4)
% hold on
% plot(XE(:,L+1),YE(:,L+1),'.','MarkerEdgeColor','r')
% % hold on
% % contour(dmx,dmy,dIS,'k');
% % hold off
% % %title('NUMTRA=2;T=50; L=9000; x0=24.1; y0=17')
% % plot(XE(1,1),YE(1,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','r')
% % for i=1:1:NUMTRA
% %     plot(XE(i,:),YE(i,:),'k');
%     plot(XE(i,L+1),YE(i,L+1),'rs','MarkerEdgeColor','k','MarkerFaceColor','g')
% end
% figure(2)
% hold on
% plot(XE(:,L+1),YE(:,L+1),'rs','MarkerEdgeColor','k','MarkerFaceColor','g')
% %title('Densidad por Nucleo Gaussiano con: NUMTRA=2;T=50; L=9000; x0=24.1; y0=17')
% plot(XE(1,1),YE(1,1),'rs','MarkerEdgeColor','k','MarkerFaceColor','r')
% contour(domx,domy,H,[0 0],'black');
% surf(xd,yd,dK);view(0,90);
%%[XEDBM,YEDBM, XEm, YEm,contDBM,L,bvnd,psm,psb]=TDBM(H,DXIS,DYIS,dIS,dISx,dISy,PIXTDBM,PIYTDBM,L); 
end