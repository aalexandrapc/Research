function [dISxPOS,dISyPOS,dISPOS]=interpolacionesDBM(dISx,dISy,dIS,ICELT,JCELT,PARTX,PARTY);       %% [dISxPOS,dISyPOS,dISPOS]=interpolacionesDBM(H,dISx,dISy,dIS,ICELT,JCELT,PARTX,PARTY);

        if (dIS(ICELT,JCELT)~=0)         
            dIS01 = (1-PARTX)*dIS(ICELT,JCELT)+PARTX*dIS(ICELT, JCELT+1);
            dIS32 = (1-PARTX)*dIS(ICELT+1,JCELT)+PARTX*dIS(ICELT+1,JCELT+1);
            dISPOS = (1-PARTY)*dIS01+PARTY*dIS32;
            
            dISx01 = (1-PARTX)*dISx(ICELT,JCELT)+PARTX*dISx(ICELT, JCELT+1);
            dISx32 = (1-PARTX)*dISx(ICELT+1,JCELT)+PARTX*dISx(ICELT+1,JCELT+1);
            dISxPOS = (1-PARTY)*dISx01+PARTY*dISx32;
            
            dISy01 = (1-PARTX)*dISy(ICELT,JCELT)+PARTX*dISy(ICELT, JCELT+1);
            dISy32 = (1-PARTX)*dISy(ICELT+1,JCELT)+PARTX*dISy(ICELT+1,JCELT+1);
            dISyPOS = (1-PARTY)*dISy01+PARTY*dISy32;
            
        else
               dISPOS = 0;
               dISxPOS = 0;
               dISyPOS = 0;
        end
end