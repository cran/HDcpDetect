        subroutine code1(n,m,h1,h2,a,c)
        implicit double precision (a-h,p-z)
        integer n,m,h1,h2,h3,h4
        integer CT

        double precision a(n,n)
        
        c=0.0
        CT=0
        if(h1.GE.0) then
        do 100 i=1,n-h1
         if (h2 .GE.0) then     
           do 200 j=1,n-h2
             if((j-i-h1)>m .OR. (i-j-h2)>m) then
             c=c+a(i,j+h2)*a(i+h1,j)
             CT=CT+1   
             else
             c=c+0.0
            end if
200    continue
         else
           h3=(-1)*h2 
           do 300 j=(h3+1),n
             if((i-j)>m .OR.(j-h3-i-h1)>m) then
             c=c+a(i,j-h3)*a(i+h1,j)
             CT=CT+1   
             else
             c=c+0.0
            end if
300    continue
         end if        
100    continue
       else
       h4=(-1)*h1
        do 400 i=(h4+1),n
         if (h2 .GE.0) then     
           do 500 j=1,n-h2
             if((j-i)>m .OR.(i-h4-j-h2)>m) then
             c=c+a(i,j+h2)*a(i-h4,j)
             CT=CT+1   
             else
             c=c+0.0
            end if
500    continue
         else
           h3=(-1)*h2 
           do 600 j=(h3+1),n
             if((j-h3-i)>m .OR.(i-h4-j)>m) then
             c=c+a(i,j-h3)*a(i-h4,j)
             CT=CT+1   
             else
             c=c+0.0
            end if
600    continue
         end if        
400    continue
       end if       
       c=c/CT        
       return 
       end
        
        
