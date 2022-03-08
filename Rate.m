function f=Rate(Reak,tBegin)

T=280;
N2=1.960e19;
O2=5.100e18;
H2O=5.100D17;
Pres=850.0;
p0=1013.25;
mAir = (N2+O2) * 298.15e0 /T * Pres / p0;

switch Reak.Type
    
    case 'PHOTO'
        suntime=mod(tBegin/3600,24);% REAL(dp), PARAMETER :: SunRise=4.50_dp, SunSet=19.50_dp
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO2'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO3'
        suntime=mod(tBegin/3600,24);
        if suntime>=4.50 && suntime<=19.50
            sun=updatesun(suntime);
            Meff=1;
            k=Reak.con(1)*sun*sun*sun;
            f=k*Meff;
        else
            f=0;
        end
    case 'PHOTO1'
        disp('type unknown')
    case 'PHOTMCM'
        %         chi(1)=Zenith(tBegin);
        %         chi(2)=1/COS(chi(1));
        %
        %
        %         IF ( chi(1) < PiHalf ) THEN
        %           ChiZmcm  = EXP( -iR%PHOTmcm(:,3) * chi(2) )
        %           yChiZmcm = chi(1) ** iR%PHOTmcm(:,2)
        %           k(iR%iPHOTmcm) = Dust * iR%PHOTmcm(:,1) * yChiZmcm * ChiZmcm
        %         END IF
        chi(1)=Zenith(tBegin);
        chi(2)=1/cos(chi(1));
        if chi(1)<pi/2
            ChiZmcm=exp(-Reak.con(3) * chi(2));
            yChiZmcm = chi(1)^Reak.con(2);
            f=Reak.con(1)*yChiZmcm*ChiZmcm;
        else
            f=0;
        end
        
    case 'SPEC3'
        %         k1s = iR%SPEC3(:,1)*EXP(iR%SPEC3(:,2)*T(6))
        %         k2s = iR%SPEC3(:,3)*EXP(iR%SPEC3(:,4)*T(6))
        %         k3s = iR%SPEC3(:,5)*EXP(iR%SPEC3(:,6)*T(6))*mAir
        %         k(iR%iSPEC3) = k1s+k3s/(One+k3s/k2s)
        k1s=Reak.con(1)*exp(Reak.con(2)/T);
        k2s=Reak.con(3)*exp(Reak.con(4)/T);
        k3s=Reak.con(5)*exp(Reak.con(6)/T)*mAir;
        f=k1s+k3s/(1+k3s/k2s);
    case 'SPEC2MCM'
        f=Reak.con(1)*(T/300)^Reak.con(2)*exp(Reak.con(3)/T);
    case 'SPEC3MCM'
        %         k(iR%iSPEC3mcm) = iR%SPEC3mcm(:,1)*(ONE+mAir/iR%SPEC3mcm(:,2))
        f=Reak.con(1)*(1+mAir/Reak.con(2));
    case 'SPEC4MCM'
        %         k(iR%iSPEC4mcm) = iR%SPEC4mcm(:,1)*(ONE+iR%SPEC4mcm(:,2) &
        %         &               * EXP(iR%SPEC4mcm(:,3)*T(6))*H2O)*EXP(iR%SPEC4mcm(:,4)*T(6))
        f=Reak.con(1)*(1+Reak.con(2)*exp(Reak.con(3)/T)*H2O)*exp(Reak.con(4)/T);
    case 'SPEC6MCM'
        k1=Reak.con(1)*exp(Reak.con(2)/T);
        k2=Reak.con(3)*exp(Reak.con(4)/T);
        f=k1*(1-k2);
    case 'SPEC7MCM'
        k1=Reak.con(1)*exp(Reak.con(2)/T);
        k2=Reak.con(3)*exp(Reak.con(4)/T);
        f=k1*(Reak.con(5)-Reak.con(6)/(1+k2));
    case 'TROEMCM'
        %         InvRefTemp=1.0e0/298.15e0;
        %         k1mcm = iR%TROEmcm(:,1)*(T(1)*InvRefTemp)**iR%TROEmcm(:,2)*EXP(iR%TROEmcm(:,3)*T(6))*mAir
        %         k2mcm = iR%TROEmcm(:,4)*(T(1)*InvRefTemp)**iR%TROEmcm(:,5)*EXP(iR%TROEmcm(:,6)*T(6))
        %         Fc    = iR%TROEmcm(:,7)*EXP(iR%TROEmcm(:,8)*T(6))+iR%TROEmcm(:,9)*EXP(T(1)/iR%TROEmcm(:,10))
        %         tmpTROE = LOG10(k1mcm/k2mcm)/(n-d*LOG10(Fc))
        %         k(iR%iTROEmcm) = k1mcm/(One+k1mcm/k2mcm)*Fc**(One/(One+tmpTROE*tmpTROE))
        %         n = 0.75_dp, d = 1.27_dp
        k1=Reak.con(1)*(T/298.15e0)^Reak.con(2)*exp(Reak.con(3)/T)*mAir;
        k2=Reak.con(4)*(T/298.15e0)^Reak.con(5)*exp(Reak.con(6)/T);
        Fc=Reak.con(7)*exp(Reak.con(8)/T)+Reak.con(9)*exp(T/Reak.con(10));
        tmp=log10(k1/k2)/(0.75-1.27*log10(Fc));
        f=k1/(1+k1/k2)*Fc^(1/(1+tmp*tmp));
    case 'CONST'
        f=Reak.con(1);
    case 'TEMP1'
        f=Reak.con(1)*exp(-Reak.con(2)/T);
    otherwise
        f=0;
        disp('type unknown')
        Reak.Type
        
end
switch Reak.Factor
    case '$O2O2'
        f=f*O2*O2;
    case '$O2N2'
        f=f*O2*N2;
    case '$O2'
        f=f*O2;
    case '$N2'
        f=f*N2;
    case '$H2O'
        f=f*H2O;
    case '$M'
        f=f*mAir;
    case '$RO2'
        f=f;
    case ''
    otherwise
       Reak.Factor 
end
end

function sun=updatesun(Tlocal)
%calculation sun for SmallStratoKPP
SunRise=4.5;
SunSet=19.5;
Ttmp = (2*Tlocal-SunRise-SunSet) / (SunSet-SunRise);
if Ttmp>0
    Ttmp=Ttmp*Ttmp;
else
    Ttmp=-Ttmp*Ttmp;
end
sun = (1+cos(pi*Ttmp)) * 0.5;
end

