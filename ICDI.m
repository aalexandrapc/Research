function [NDVIT,BwT,DXIS,DYIS,cbH,fbH,cb,fb]= ICDI(H); 
%Esta función calcula el índice de clorofila NDVI de las imagenes satelitales correspondientes a una fecha
%pero la que sale es NDVIT que es la matriz cuyas filas  son las filas de NDVI ordenadas desde la última hasta
%la primenra
%B6=imread('C:\Users\Angie\Documents\AngiePineda\DoctoradoAngie\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B6.TIF', ...
     %'PixelRegion', {[3294 4000], [683 1915]});  [filas,columnas] con esta
     %instruccion selee y se corta la imagen satelital
     %[3294 4000], [740 1890] este recorte funciona pero puede traer
     %problemas con lossaltos del Brawniano .....
     
B6=imread('C:\Users\Angie\Desktop\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B6.TIF', ...
'PixelRegion', {[3305 3980], [750 1880]});                                 %[3300 3983], [750 1885] [3294 4000], [750 1885] 740 1890 con esta instruccion selee y se corta la imagen satelital/ {[filas],[columnas]} /{[3294 4000],[735 1915]}                              %
figure (30) ,imshow(B6)                                                                           
Bw = B6 < 9155;
figure (1), imshow(Bw), title('Identificación Lago banda B6');
Bw =  bwselect(Bw, 424, 280, 4);
figure (2), imshow(Bw), title('Mascara en función de B6');

B4=imread('C:\Users\Angie\Desktop\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B4.TIF', ...
     'PixelRegion', {[3305 3980], [750 1880]});  
B5=imread('C:\Users\Angie\Desktop\IMAGENES PARA COEFICIENTE DE DIFUSIÓN\2014-05-13T171559Z\LC80040532014133LGN00_B5.TIF', ...
     'PixelRegion', {[3305 3980], [750 1880]});
figure (3), imshow(B4,[]), title('Identificación Lago banda  B4');
figure (4), imshow(B5,[]), title('Identificación  Lago banda  B5');
B4(~Bw)=0;                                                                 % se hacen 0 los valores de la matriz de B4 que no esten en la zona limitada por  la mascara
B5(~Bw)=0;                                                                 % se hacen 0 los valores de la matriz de B5 que no esten en la zona limitada por  la mascara
figure(5), imshow(B4,[]),title('INDICE B4');
figure(6), imshow(B4,[]),title('INDICE B5');
%[f4,c4]=size(B4);
%%Ajuste de malla
% %% Calculo de dimensiones del lago en función de H
% H(~(H~=4.5))=0;                                                            % se hacen 0 las entradas de H correspondientes a profundidad 0(en el codigo aparece 4.5) en malla "batimetria_valencia_Saa"
% cbH=find(any(H));                                                          % indices de las columnas de H que tienen entradas distintasde 0
% fbH=find(any(H'));                                                         % indices de las filas de H que tienen entradas distintas a 0
% LLV= (1/3)*(length(cbH));                                                  % largo del lagode Valencia en función de malla "batimetria_valencia_Saa"
% ALV= (1/3)*(length(fbH));                                                  % Ancho  del lagode Valencia en función de malla "batimetria_valencia_Saa"
% %% Cálculo de los parámetros de la malla rectangular correspondiente a la imagen satelital
% cb=find(any(B4));                                                          % indices de las columnas de banda B4 que tienen entradas distintasde 0
% fb=find(any(B4'));                                                         % indices de las filas debanda B4  que tienen entradas distintasde 0
% DXIS=LLV/(length(cb));                                                    % longitud horizontal de un rectangulo correspondiente a la malla de laimagen satelital
% DYIS=ALV/(length(fb));                                                    % longitud vertical de un rectangulo correspondiente a la malla de laimagen satelital
% UU=any(B4);
% VV=any(B4');
%% Cálculo del índice de vegetación
 
NDVI = zeros(size(B4),'single');

RM4 = 2.0000E-05;
RM5 = 2.0000E-05;
RA4 = -0.100000;
RA5 = -0.100000;

RL4=RM4*single(B4(Bw))+RA4;
RL5=RM5*single(B5(Bw))+RA5;

NDVI(Bw) = (single(B5(Bw))-single(B4(Bw)))./(single(B5(Bw))+single(B4(Bw))+0.1);
NDVI(Bw) = (single(RL5)-single(RL4))./(single(RL5)+single(RL4)+0.1);
figure (7), imshow(NDVI,[]), title('Indice de vegetación en el Lago de Valencia')
[nrN,ncN]=size(NDVI);                                                      %Numero de filas y columnas  de NDVI        
NDVIT=zeros(size(NDVI));
BwT=zeros(size(Bw));
for i=1:nrN
    NDVIT(i,:)=NDVI(nrN+1-i,:);                                            %NDVIT es una matriz que ordena las filas de NDVI desde la ultima hasta la primera
    BwT(i,:)=Bw(nrN+1-i,:);
end
figure (8), imshow(NDVIT,[]), title('Indice de Clorofila en el lago volteado') 
%% Calculo de dimensiones del lago en función de H
H(~(H~=4.5))=0;                                                            % se hacen 0 las entradas de H correspondientes a profundidad 0(en el codigo aparece 4.5) en malla "batimetria_valencia_Saa"
cbH=find(any(H));                                                          % indices de las columnas de H que tienen entradas distintasde 0
fbH=find(any(H'));                                                         % indices de las filas de H que tienen entradas distintas a 0
LLV= (1/3)*(length(cbH));                                                  % largo del lagode Valencia en función de malla "batimetria_valencia_Saa"
ALV= (1/3)*(length(fbH));                                                  % Ancho  del lagode Valencia en función de malla "batimetria_valencia_Saa"
%% Cálculo de los parámetros de la malla rectangular correspondiente a la imagen satelital
cb=find(any(BwT(2:end-1,2:end-1)));                                                         % indices de las columnas de banda B4 que tienen entradas distintasde 0
fb=find(any((BwT(2:end-1,2:end-1))'));                                                        % indices de las filas debanda B4  que tienen entradas distintasde 0
DXIS=LLV/(length(cb));%DXIS=LLV/(length(cb));                              % longitud horizontal de un rectangulo correspondiente a la malla de laimagen satelital
DYIS=ALV/(length(fb)); %DYIS=ALV/(length(fb));                             % longitud vertical de un rectangulo correspondiente a la malla de laimagen satelital
% UU=any(BwT);
% VV=any(BwT');
end