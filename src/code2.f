        subroutine code2(n,m,h,a,c)
        implicit double precision (a-h,p-z)
        integer n,m,h
        integer CT
        
        double precision a(n,n)
        
        c=0.0
        CT=0 
        do 100 i=1,n 
         do 200 j=1,n-h
          if ((j-i)>m .OR.(i-j-h)>m) then        
            do 300 k=1,n   
         if (ABS(k-i)>m.AND.ABS(k-j)>m.AND.ABS(k-j-h)>m) then
                 c=c+a(i,j)*a(j+h,k)
                 CT=CT+1
               end if
                
300    continue
            end if
200    continue        
100    continue        
       c=c/CT
       return 
       end
       
       
