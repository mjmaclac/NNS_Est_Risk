out_num=150;

t_future=repmat([1:(len+out_num)]',1,len+out_num); u_future=repmat([1:(len+out_num)],len+out_num,1);

find_val=find(fval_array_temp==min(fval_array_temp(:)));

[I1,I2,I3]=ind2sub([profs_num profs_num profs_num],find_val);

a0_imp=[a_0(I1) a_1(I2) a_2(I3)];

if a0(2)~=0
    pr_future=tril(1 ./ (1+exp(a0_imp(1)-a0_imp(2)*(t_future-u_future).^2)));
else
    pr_future=tril(1 ./ (1+exp(a0(1))));
end

pr_intro_u_NOTdisc_b4t = tril(cumprod(1 - pr_future ,1)); 
y=pr_intro_u_NOTdisc_b4t;

yy=zeros(size(y,1),1);
yy_10=zeros(size(y,1),1);
yy_90=zeros(size(y,1),1);
for i=1:size(y,1)
    for j=1:size(y,2)
        if y(i,j)<=1 && y(i,j)>=0.5 && y(i,j)~=0
            yy(j,1)=yy(j,1)+1; 
        end
        if y(i,j)<=1 && y(i,j)>=0.025 && y(i,j)~=0
            yy_10(j,1)=yy_10(j,1)+1; 
        end
        if y(i,j)<=1 && y(i,j)>=0.975 && y(i,j)~=0
            yy_90(j,1)=yy_90(j,1)+1; 
        end
    end
end

[median(yy_90) median(yy) median(yy_10)]
