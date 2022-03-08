  function Chi=Zenith(Time)
    %
    %----------------------------------------------------------------------!
    %
    IMN=[31,28,31,30,31,30,31,31,30,31,30,31];
    ONE=1;
    TWO=2;
    rTWO=ONE/TWO;
    THREE=3;
    FOUR=4;
    rFOUR=ONE/FOUR;
    FIVE=5;
    HOUR=3600.0;
    rlat=45.00; %50.65;
    rlon=0;     %10.77;
    DR=pi/180.0;
    iDate=010621;
    % set GMT
    GMT = Time / HOUR;
    %
    %  convert to radians
    RLT = rlat*DR;
    RPHI = rlon*DR;
    %
    %  parse date
    IIYEAR = round(iDate/10000);
    IYEAR = 19*100 + IIYEAR;
    if IIYEAR <= 50 
        IYEAR = IYEAR + 100; 
    end    
    IMTH = round((iDate - IIYEAR*10000)/100);
    IDAY = iDate - IIYEAR*10000 - IMTH*100;
    %
    %  identify and correct leap years
    IIY = round((IIYEAR/4))*4;
    if IIY==IIYEAR 
        IMN(2) = 29;
    end     
    %
    %  count days from Dec.31,1973 to Jan 1, YEAR, then add to 2,442,047.5
    YREF =  2442047.5e0;
    NYEARS = IYEAR - 1974;
    LEAP = round((NYEARS+1)/4);
    if NYEARS<=-1 
       LEAP = round((NYEARS-2)/4);
    end   
    NOLEAP = NYEARS - LEAP;
    YR = YREF + 365.0e0*NOLEAP + 366.0e0*LEAP;
    %
    IJD = 0;
    IN = IMTH - 1;
    if IN~=0 
        for I=1:IN
            IJD=IJD+IMN(I);
        end
    IJD = IJD + IDAY;
    else
     IJD = IDAY;
   end  
     IJ = IYEAR - 1973;


%   IF(IN.EQ.0) GO TO 40
%   DO 30 I=1,IN
%   IJD = IJD + IMN(I)
% 30   CONTINUE
%   IJD = IJD + IDAY
%   GO TO 50
% 40   IJD = IDAY
% 50   IJ = IYEAR - 1973
    %
    %      print julian days current "ijd"
    JD = IJD + (YR - YREF);
    D = JD + GMT/24.0e0;
    %
    %      calc geom mean longitude
    ML = 279.2801988e0 + .9856473354e0*D + 2.267e-13*D*D;
    RML = ML*DR;
    %
    %      calc equation of time in sec
    %      w = mean long of perigee
    %      e = eccentricity
    %      epsi = mean obliquity of ecliptic
    W = 282.4932328e0 + 4.70684e-5*D + 3.39e-13*D*D;
    WR = W*DR;
    EC = 1.6720041e-2 - 1.1444e-9*D - 9.4e-17*D*D;
    EPSI = 23.44266511 - 3.5626e-7*D - 1.23e-15*D*D;
    PEPSI = EPSI*DR;
    YT = (tan(PEPSI*rTWO))^2;
    CW = cos(WR);
    SW = sin(WR);
    SSW = sin(TWO*WR);
    EYT = TWO*EC*YT;
    FEQT1 = sin(RML)*(-EYT*CW - TWO*EC*CW);
    FEQT2 = cos(RML)*(TWO*EC*SW - EYT*SW);
    FEQT3 = sin(TWO*RML)*(YT - (FIVE*EC*EC*rFOUR)*(CW*CW-SW*SW));
    FEQT4 = cos(TWO*RML)*(FIVE*EC^2*SSW*rFOUR);
    FEQT5 = sin(THREE*RML)*(EYT*CW);
    FEQT6 = cos(THREE*RML)*(-EYT*SW);
    FEQT7 = -sin(FOUR*RML)*(rTWO*YT*YT);
    FEQT = FEQT1 + FEQT2 + FEQT3 + FEQT4 + FEQT5 + FEQT6 + FEQT7;
    EQT = FEQT*13751.0e0;
    %
    %   convert eq of time from sec to deg
    REQT = EQT/240.0e0;
    %
    %   calc right ascension in rads
    RA = ML - REQT;
    RRA = RA*DR;
    %
    %   calc declination in rads, deg
    TAB = 0.43360e0*sin(RRA);
    RDECL = atan(TAB);
    DECL = RDECL/DR;
    %
    %   calc local hour angle
    LBGMT = 12.0e0 - EQT/HOUR + rlon*24.0e0/360.0e0;
    LZGMT = 15.0e0*(GMT - LBGMT);
    ZPT = LZGMT*DR;
    CSZ = sin(RLT)*sin(RDECL) + cos(RLT)*cos(RDECL)*cos(ZPT);
    ZR = acos(CSZ);
    % 
    %   calc local solar azimuth
    CAZ = (sin(RDECL) - sin(RLT)*cos(ZR))/(cos(RLT)*sin(ZR));
    RAZ = acos(CAZ);
    AZIMUTH = RAZ/DR;
    %
    %--- set Zenith Angle
    Chi =  1.745329252e-02 * ZR/DR;
    end
