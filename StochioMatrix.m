function [A] = StochioMatrix(Reak,Species)

A=zeros(size(Reak,1),size(Species,1));
for l=1:size(Reak,1)
    for j=1:Reak(l).Right
        PosR=GetPos(Reak(l).NameR(j).Name,Species);
        A(l,PosR)=A(l,PosR)+Reak(l).KoeffR(j);
    end
    for j=1:Reak(l).Left
        PosL=GetPos(Reak(l).NameL(j).Name,Species);
        A(l,PosL)=A(l,PosL)-Reak(l).KoeffL(j);
    end
end
A=sparse(A);
end

