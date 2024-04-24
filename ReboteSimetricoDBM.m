function [XE,YE]=ReboteSimetricoDBM(i,j,dIS,DX,DY,XPOS,YPOS,TRACER,XE,YE);  
       %%[XE,YE,psb]=controldebordeDBM(i,j,dIS,ICELT,JCELT,DX,DY,XPOS,YPOS,TRACER,XE,YE,psb);
         JCELP = length((DX:DX:XPOS));                                     %% índice j  correspondiente al punto (XPOS,YPOS), en la matriz  determinada por la malla de las velocidades del fluido                                   
         ICELP = length((DY:DY:YPOS));                                     %% índice i  correspondiente al punto (XPOS,YPOS), en la matriz  determinada por la malla de las velocidades del fluido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dIS(ICELP,JCELP)~=0)                                                    %% Condición que verifican los puntos sobre la superficie del lago
    XE(i,j+1)=XPOS;
    YE(i,j+1)=YPOS;
else                                                                       %% cálculo del rebote simétrico para particulas que salen del dominio del lago
    dpb=zeros(1,4);                                                        %% Esta matriz será llenada con las distancias del punto al borde
    l=zeros(1,4);                                                          %% Sus entradas registrarán el número de celdas con agua, a la izquierda, derecha, abajo y arriba del punto fuera del dominio
    C1=find(dIS(ICELP,1:JCELP)~=0);                                        %% indices j de las celdas  con agua que estan a la izquierda del punto fuera del dominio
    l(1)=length(C1);                                                       %% número de celdas a la izquierda con agua
    dIS1=dIS;
    dIS1(ICELP,1:JCELP)=zeros(1,JCELP);
    C2=find(dIS1(ICELP,:)~=0);                                             %% celdas a la derecha con agua
    l(2)=length(C2);                                                       %% número de celdas a la derecha con agua
    F1=find(dIS(1:ICELP,JCELP)~=0);                                        %% celdas  hacia abajo con agua 
    l(3)=length(F1);                                                       %% número de filas hacia abajo con agua
    dIS2=dIS;
    dIS2(1:ICELP,JCELP)=zeros(ICELP,1);
    F2=find(dIS2(:,JCELP)~=0);                                             %% celdas hacia arriba con agua
    l(4)=length(F2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (l(1)>=1)
        dpb(1)=abs(XPOS-(C1(end)+1)*DX);                                   
    else
        dpb(1)=3000;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (l(2)>=1)
        dpb(2)=abs((C2(1)*DX)-XPOS);
    else
        dpb(2)=3000;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (l(3)>=1)
        dpb(3)=abs(YPOS-(F1(end)+1)*DY);
    else
        dpb(3)=3000;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (l(4)>=1)
        dpb(4)=abs((F2(1)*DY)-YPOS);
    else
        dpb(4)=3000;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mdpb=min(dpb);                                                         %% mínimo de las distancias del punto hacia los bordes del lago que están a su alrededor
    if (mdpb<=2999)                                                       
        Dis=find(dpb == mdpb);                                             %% índices de las distancias mínimas al borde del lago                                                                                     
        lD=length(Dis);
        XPOSp=zeros(1,lD);
        YPOSp=zeros(1,lD);
        ICELPp=zeros(1,lD);
        JCELPp=zeros(1,lD);
        dISp=zeros(1,lD);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cálculo de los posibles valores de XE(i,j+1) e Y(i,j+1)
        for h=1:lD
            if (Dis(h)== 1)
               XPOSp(h)=((C1(end)+1)*DX)-dpb(1);
               YPOSp(h)=YPOS;
               ICELPp(h)= ICELP;
               JCELPp(h)= length((DX:DX:XPOSp(h)));
               dISp(h)=dIS(ICELPp(h),JCELPp(h));
            elseif (Dis(h)== 2)
                XPOSp(h)= (C2(1)*DX)+dpb(2);
                YPOSp(h)= YPOS;
                ICELPp(h)= ICELP;
                JCELPp(h)= length((DX:DX:XPOSp(h)));
                dISp(h)=dIS(ICELPp(h),JCELPp(h));
            elseif (Dis(h)== 3)
                XPOSp(h)=XPOS;
                YPOSp(h)=((F1(end)+1)*DY)-dpb(3);
                ICELPp(h)=length((DY:DY:YPOSp(h))) ;
                JCELPp(h)=JCELP;
                dISp(h)=dIS(ICELPp(h),JCELPp(h));
            else
                XPOSp(h)=XPOS;
                YPOSp(h)=(F2(1)*DY)+dpb(4);
                ICELPp(h)=length((DY:DY:YPOSp(h)));
                JCELPp(h)=JCELP;
                dISp(h)=dIS(ICELPp(h),JCELPp(h));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CELF=find(dISp ~= 0);                                                %% Con esta instrucción se verifica si, luego del rebote simétrico, el punto está en una celda con agua 
        lCELF=length(CELF);
        if (lCELF >= 1)
            IndiceFi=CELF(1);
            XE(i,j+1)= XPOSp(IndiceFi);
            YE(i,j+1)= YPOSp(IndiceFi);
        else
           XE(i,j+1) = TRACER(1);
           YE(i,j+1) = TRACER(2);
        end
    else
        XE(i,j+1)= TRACER(1);
        YE(i,j+1)= TRACER(2);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
