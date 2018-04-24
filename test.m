User=[1,2,3,4,5];
beta=2.5;
User_SBSid=zeros(1,5);
a=1;
for i=1:5
   if User(i)>beta 
      User_SBSid(a)=User(i);
      a=a+1;
   end
end