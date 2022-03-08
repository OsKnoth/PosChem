function M=Jacobian(t,y)
global Reak Species
M=zeros(size(Species,1),size(Species,1));
for l=1:size(Reak,1)
    r=Rate(Reak(l),t);
    for i=1:Reak(l).Left
        PosL=GetPos(Reak(l).NameL(i).Name,Species);
        rr=r;
        for j=1:Reak(l).Left
            if i~=j
                Pos=GetPos(Reak(l).NameL(j).Name,Species);
                rr=rr*y(Pos);
            end
        end
        for j=1:Reak(l).Left
            PosR=GetPos(Reak(l).NameL(j).Name,Species);
            M(PosR,PosL)=M(PosR,PosL)-Reak(l).KoeffL(j)*rr;
        end
        for j=1:Reak(l).Right
            PosR=GetPos(Reak(l).NameR(j).Name,Species);
            M(PosR,PosL)=M(PosR,PosL)+Reak(l).KoeffR(j)*rr;
        end
    end
end
M=sparse(M);
end


