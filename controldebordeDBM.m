function [XE,YE,psb]=controldebordeDBM(i,j,dIS,ICELT,JCELT,DX,DY,XPOS,YPOS,TRACER,XE,YE,psb);
        JCELP = length((DX:DX:XPOS));
        ICELP = length((DY:DY:YPOS));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (dIS(ICELP,JCELP)>0)                                                                                      
    XE(i,j+1)=XPOS; 
    YE(i,j+1)=YPOS;
else
    psb=psb+1;                                                             %%psb cuenta el númerode veces que una particula rebota desde afuera hacia adentro del lago
    DISTX=XPOS-TRACER(1);                                                  
    DISTY=YPOS-TRACER(2);
    if ((DISTX == 0) && (DISTY == 0))
        ALPHA = 0;
    elseif ((DISTX == 0) && (DISTY ~= 0))    
           if (DISTY > 0) 
                ALPHA = (pi/2);
           else
             ALPHA = (3*pi/2);
           end 
    elseif ((DISTY == 0) && (DISTX ~= 0))
            if (DISTX > 0)
            ALPHA = 0;
            else
            ALPHA = pi; 
            end
    else 
          ALPHA=atan2(DISTY,DISTX);
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ALPHA < 0)
        ALPHA=ALPHA+(2*pi);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if (ALPHA == 0)
        epsilon=((JCELT+1)*DX-TRACER(1))/3;
        XE(i,j+1)= (JCELT+1)*DX-epsilon;                                                        
        YE(i,j+1)=YPOS;
    end
    if ((ALPHA > 0) && (ALPHA < (pi/2)))
         Y=((TRACER(2)-YPOS)/(TRACER(1)-XPOS))*((JCELT+1)*DX-XPOS)+YPOS;
         X=((TRACER(1)-XPOS)/(TRACER(2)-YPOS))*((ICELT+1)*DY-YPOS)+XPOS;
         if (Y < (ICELT+1)*DY)
             epsilon=((JCELT+1)*DX-TRACER(1))/3; 
             dist=tan(ALPHA)*epsilon;
             XE(i,j+1)=((JCELT+1)*DX)-epsilon;
             YE(i,j+1)=Y-dist; 
         elseif (X < (JCELT+1)*DX)
                epsilon=((ICELT+1)*DY-TRACER(2))/3; 
                dist=epsilon/tan(ALPHA);
                XE(i,j+1)=X-dist;
                YE(i,j+1)=(ICELT +1)*DY-epsilon;
         else                                                       
                epsilon=((ICELT+1)*DY-TRACER(2))/3; 
                dist=epsilon/tan(ALPHA);
                XE(i,j+1)=(JCELT+1)*DX-dist;
                YE(i,j+1)=(ICELT+1)*DY-epsilon;
         end
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    if (ALPHA == (pi/2))
        epsilon=((ICELT+1)*DY-TRACER(2))/3;
        XE(i,j+1)= XPOS;
        YE(i,j+1)= (ICELT+1)*DY-epsilon;
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (((pi/2) < ALPHA ) && (ALPHA < pi))
                THETA = pi-ALPHA;
                Y=((TRACER(2)-YPOS)/(TRACER(1)-XPOS))*(JCELT*DX -XPOS)+YPOS;
                X=((TRACER(1)-XPOS)/(TRACER(2)-YPOS))*((ICELT+1)*DY-YPOS)+XPOS;
                if (Y < (ICELT+1)*DY)
                    epsilon=(TRACER(1)-JCELT*DX)/3; 
                    dist=tan(THETA)*epsilon;
                    XE(i,j+1)=JCELT*DX+epsilon;
                    YE(i,j+1)=Y-dist;
                elseif (X > JCELT*DX)
                    epsilon=((ICELT+1)*DY-TRACER(2))/3; 
                    dist=epsilon/tan(THETA);
                    XE(i,j+1)=X+dist;
                    YE(i,j+1)=(ICELT+1)*DY-epsilon;
                else                                                       %justo aqui tengo la duda de colocar o no else o elseif junto con la condición correspondiente Y=ICELT+1 
                    epsilon=(TRACER(1)-JCELT*DX)/3; 
                    dist=epsilon*tan(THETA);
                    XE(i,j+1)=JCELT*DX+epsilon;
                    YE(i,j+1)=(ICELT+1)*DY-dist;
                end  
            end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
            if (ALPHA == pi)
                epsilon=(TRACER(1)- JCELT*DX)/3;
                XE(i,j+1)=JCELT*DX + epsilon;
                YE(i,j+1)=YPOS;
            end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            if ((pi < ALPHA) && (ALPHA < (3*pi/2)))
                THETA = ALPHA - pi;
                Y=((TRACER(2)-YPOS)/(TRACER(1)-XPOS))*(JCELT*DX-XPOS)+YPOS;
                X=((TRACER(1)-XPOS)/(TRACER(2)-YPOS))*(ICELT*DY-YPOS)+XPOS;
                if (Y > ICELT*DY)
                    epsilon=(TRACER(1)-JCELT*DX)/3; 
                    dist=tan(THETA)*epsilon;
                    XE(i,j+1)=JCELT*DX +epsilon;
                    YE(i,j+1)=Y+dist;
                elseif (X > JCELT*DX)
                    epsilon=(TRACER(2)-ICELT*DY)/3; 
                    dist=epsilon/tan(THETA);
                    XE(i,j+1)=X+dist;
                    YE(i,j+1)=ICELT*DY+epsilon;
                else                                                        
                    epsilon=(TRACER(1)-JCELT*DX)/3; 
                    dist=epsilon*tan(THETA);
                    XE(i,j+1)=JCELT*DX +epsilon;
                    YE(i,j+1)=ICELT*DY +dist;
                end
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if (ALPHA == (3*pi/2))
                epsilon=(TRACER(2)-ICELT*DY)/3;
                XE(i,j+1)=XPOS;
                YE(i,j+1)=ICELT*DY+epsilon;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            if (((3*pi/2) < ALPHA) && (ALPHA < (2*pi)))
                THETA = (2*pi-ALPHA);
                Y=((TRACER(2)-YPOS)/(TRACER(1)-XPOS))*((JCELT+1)*DX -XPOS)+YPOS;
                X=((TRACER(1)-XPOS)/(TRACER(2)-YPOS))*(ICELT*DY-YPOS)+XPOS;
                if (Y > ICELT*DY)
                    epsilon=((JCELT + 1)*DX-TRACER(1))/3; 
                    dist=tan(THETA)*epsilon;
                    XE(i,j+1)=(JCELT+1)*DX-epsilon;
                    YE(i,j+1)=Y+dist;
                elseif (X < (JCELT+1)*DX )
                    epsilon=(TRACER(2)-ICELT*DY)/3; 
                    dist=epsilon/tan(THETA);
                    XE(i,j+1)=X-dist;
                    YE(i,j+1)=ICELT*DY + epsilon;
                else                                                       
                    epsilon=((JCELT+1)*DX -TRACER(1))/3; 
                    dist=epsilon*tan(THETA);
                    XE(i,j+1)=(JCELT+1)*DX-epsilon;
                    YE(i,j+1)=ICELT*DY+dist;
                end
            end
 end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
        