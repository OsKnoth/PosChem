function M=MHatMatrix(Reak,Species,tBegin,y)
M=zeros(size(Species,1),size(Species,1));
for i=1:size(Reak,1)
  M=InsertReak(y,Reak(i),M,Species,tBegin);
end
end

