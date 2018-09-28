        subroutine code3(n,m,a,c)
        implicit double precision (a-h,p-z)
        integer n,m
        integer CT

        double precision a(n,n)
        
        c=0.0
        CT=0
        do 100 i=1,n 
         do 200 j=1,n
           if (ABS(j-i)>m) then
            do 300 k=1,n
             if (ABS(k-i)>m.AND.ABS(k-j)>m) then
              do 400 l=1,n
        if(ABS(l-i)>m.AND.ABS(l-j)>m.AND.ABS(l-k)>m) then
           c=c+a(i,j)*a(k,l)
           CT=CT+1
        end if
400    continue
             end if
300    continue
           end if
200    continue        
100    continue        
       c=c/CT
       return 
       end
