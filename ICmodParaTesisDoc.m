function [NDVIT,DXIS,DYIS,nrN,ncN,BwT,Bww1,Bww2]= ICmodParaTesisDoc
% Esta funciòn es una modificaciòn de de la funciòn IC contenida en la
% misma carpeta para revisdar el procedimiento y describirlo en la tesis doctoral
%%%%%%%%%%%%%%%%%%%%%%%%%Open mod for chapter IV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('Datacn1.mat');                                           % this files contains the initial dataset of the bathymetry H                  
[H,Hx,Hy,domx,domy]=batimetria_valencia_Saa(Datacn1);
B6=imread('C:\Users\Usuario\Desktop\CodigosLagoParaChapterIV\DistortedBrownianMotion\LC80040532014133LGN00_B6.TIF');  %lectura de B6 SIN RECORTAR                                        %con esta instruccion se lee la imagen y se recorta
figure(1), imshow(B6,[]); %, title('Identificación Lago banda  B6 sin recortar')
B6=imread('C:\Users\Usuario\Desktop\CodigosLagoParaChapterIV\DistortedBrownianMotion\LC80040532014133LGN00_B6.TIF', ... % Lectura de B6 recortada
 'PixelRegion', {[3294 4000], [740 1890]});%con esta instruccion se lee la imagen y se recorta
figure(2), imshow(B6,[]); % , title('Identificación Lago banda  B6 sin recortar')                              %
%%%%%%%%%%%%%%%%%%%%%%%% clouse mod for chapter IV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                          
Bw = B6 < 9155;
Bww1=Bw;
figure(3), imshow(Bw) %, title('Identificación Lago banda B6')
Bw =  bwselect(Bw, 424, 280, 4);
Bww2=Bw;
figure(4), imshow(Bw); % title('Mascara en función de B6')
%%%%%%%%%%%%%%%%%%%%%%%%%Open mod for chapter IV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
B4=imread('C:\Users\Usuario\Desktop\CodigosLagoParaChapterIV\DistortedBrownianMotion\LC80040532014133LGN00_B4.TIF', ...
     'PixelRegion', {[3294 4000], [740 1890]}); 
% B4=imread('C:\Users\Usuario\Desktop\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B4.TIF', ...
%      'PixelRegion', {[3294 4000], [740 1890]});  
 %%%%%%%%%%%%%%%%%%%%%%%%clouse mod for chapterIV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%Open mod for chapter IV%%%%%%%%%%%%%%%%%%%%%%%%%
B5=imread('C:\Users\Usuario\Desktop\CodigosLagoParaChapterIV\DistortedBrownianMotion\LC80040532014133LGN00_B5.TIF', ...
     'PixelRegion', {[3294 4000], [740 1890]});
%  B5=imread('C:\Users\Usuario\Desktop\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B5.TIF', ...
%      'PixelRegion', {[3294 4000], [740 1890]});
 %%%%%%%%%%%%%%%%%%%%%%%%clouse mod for chapterIV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5), imshow(B4,[]); % , title('Identificación Lago banda  B4')
figure(6), imshow(B5,[]);   %, title('Identificación  Lago banda  B5')
B4(~Bw)=0;                                                                 % se hacen 0 los valores de la matriz de B4 que no esten en la zona limitada por  la mascara
B5(~Bw)=0;                                                                 % se hacen 0 los valores de la matriz de B5 que no esten en la zona limitada por  la mascara
%figure(5), imshow(B4,[]),title('INDICE B4');
%figure(6), imshow(B4,[]),title('INDICE B5');
%[f4,c4]=size(B4);
%%Ajuste de malla
%% Calculo de dimensiones del lago en función de H
H(~(H~=4.5))=0;                                                            % se hacen 0 las entradas de H correspondientes a profundidad 0 en malla "batimetria_valencia_Saa"
cbH=find(any(H));                                                          % indices de las columnas de H que tienen entradas distintasde 0
fbH=find(any(H'));                                                         % indices de las filas de H que tienen entradas distintas a 0
cminH=cbH(1);                                                              % indice de la primera columna de H con entradas  distinto de 0                                               
fminH=fbH(1) ;                                                             % indice de la primera fila de H con entradas distintas a 0
cmaxH=cbH(end);                                                            % indice de la última  columna de H con entradas  distinto de 0
fmaxH=fbH(end);                                                            % indice de la última fila de H con entradas  distinto de 0
LLV= (1/3)*(cmaxH-cminH+1);                                                % largo del lagode Valencia en función de malla "batimetria_valencia_Saa"
ALV= (1/3)*(fmaxH-fminH+1)  ;                                              % Ancho  del lagode Valencia en función de malla "batimetria_valencia_Saa"
%% Cálculo de los parámetros de la malla rectangular correspondiente a la imagen satelital
cb=find(any(B4));                                                          % indices de las columnas de banda B4 que tienen entradas distintasde 0
fb=find(any(B4'));                                                         % indices de las filas debanda B4  que tienen entradas distintasde 0
cmin=cb(1) ;                                                               % primera columna de banda B4 con un elemento distinto de cero
fmin=fb(1) ;                                                               % primera fila de banda B4 con un elemento distinto de cero
cmax=cb(end);                                                              % última  columna de banda B4 con un elemento distinto de cero
fmax=fb(end);                                                              % última fila de banda B4 con un elemento distinto de cero
DXIS=LLV/(cmax-cmin+1);                                                    % longitud horizontal de un rectangulo correspondiente a la malla de la imagen satelital
DYIS=ALV/(fmax-fmin+1);                                                    % longitud vertical de un rectangulo correspondiente a la malla de laimagen satelital

%% Cálculo del índice de vegetación
 
NDVI = zeros(size(B4),'single');

RM4 = 2.0000E-05;                                                          % Reflactance
RM5 = 2.0000E-05;
RA4 = -0.100000;
RA5 = -0.100000;

RL4=RM4*single(B4(Bw))+RA4;
RL5=RM5*single(B5(Bw))+RA5;

NDVI(Bw) = (single(B5(Bw))-single(B4(Bw)))./(single(B5(Bw))+single(B4(Bw))+0.1);
figure(7), imshow(NDVI,[]); %, title('Indice de Clorofila en el lago volteado')
NDVI(Bw) = (single(RL5)-single(RL4))./(single(RL5)+single(RL4)+0.1);
%figure (7), imshow(NDVI,[]), title('Indice de Clorofila en el lago')
nrN=length(NDVI(:,1));                                                     %Numero de filas de NDVI
ncN=length(NDVI(1,:));                                                     %Numero de columnas de NDVI
for i=1:nrN
    NDVIT(i,:)=NDVI(nrN+1-i,:);                                            %NDVIT es una matriz que ordena las filas de NDVI desde la ultima hasta la primera
    BwT(i,:)=Bw(nrN+1-i,:);
end
figure(8), imshow(NDVIT,[]); %, title('Indice de Clorofila en el lago volteado') 
%reordenar las filas de la matriz NDVI en NDVIT permite trabajar sobre un
%lago posicionado de igual forma que en "batimetria_valencia_Saa"
end