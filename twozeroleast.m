% 
function start=twozeroleast(ss)

res=find(ss==0);
% 
if length(res)<=1
       % 
            start=100;
            
end
if length(res)==2
    if res(1)==res(2)-1
        start=res(1)-1;
       
        return;
    end
end

if length(res)>2
    k=0;
samenum=0;
for i=1:length(res)
    %
    k=i+1;
    if res(k)==res(i)+1
        samenum=samenum+1;
        if samenum>=2
           start=res(i)-samenum;
        
           break;
        else
            %
            start=100;
         
        end
    end
end
end
% 
end
