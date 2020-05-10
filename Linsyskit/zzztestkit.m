clc
disp('   ');
disp('...... Welcome! This function is to test all the m-functions in Linear Systems Toolkit ......');
disp('   ');
askfirst=input('...... Have you entered the system data (A,B,C,D) in the Workplace yet? ( 1 = yes & 0 = No )? ');
if askfirst==0
    disp('   ')
    disp('...... The computer will generate the system data with random numbers for you ......');
    disp('   ')
    sysdim=input('...... Enter a system state, input and output dimensions [n m p ] = ');
    n=sysdim(1);
    m=sysdim(2);
    p=sysdim(3);
    A=rand(n,n);
    B=rand(n,m);
    C=rand(p,n);
    disp('   ')
    disz=input('...... Choose Matrix D to be zero? 0 = Yes and 1 = No ( 0 or 1)? ');
    if disz==0
        D=C*B*0;
    else
        D=rand(p,m);
    end
    D2=D;
    C2=C;
    C1=rand(size(B,2),size(B,1));
end
A,B,C,D
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SSD ...... Hit any key to continue ......')
pause
[DD,T,nn,no,np] = ssd(A)
err=norm(DD-inv(T)*A*T)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DSSD ...... Hit any key to continue ......')
pause
[DD,T,nn,no,np] = dssd(A)
err=norm(DD-inv(T)*A*T)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function JCF ...... Hit any key to continue ......')
pause
[J,T] = jcf(A)
err=norm(J-inv(T)*A*T)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function RJD ...... Hit any key to continue ......')
pause
[J,T] = rjd(A)
err=norm(J-inv(T)*A*T)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function OSD ...... Hit any key to continue ......')
pause
[At,Ct,Ts,To,uom,Oidx] = osd(A,C)
err=norm(At-inv(Ts)*A*Ts)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function OBVIDX ...... Hit any key to continue ......')
pause
Oidx = obvidx(A,C)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function BDOSD ...... Hit any key to continue ......')
pause
[At,Ct,Ts,To,ks] = bdosd(A,C)
err=norm(At-inv(Ts)*A*Ts)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function CSD ...... Hit any key to continue ......')
pause
[At,Bt,Ts,Ti,ucm,Cidx] = csd(A,B)
err=norm(At-inv(Ts)*A*Ts)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function CTRIDX ...... Hit any key to continue ......')
pause
Cidx = ctridx(A,B)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function BDCSD ...... Hit any key to continue ......')
pause
[At,Bt,Ts,Ti,ks] = bdcsd(A,B)
err=norm(At-inv(Ts)*A*Ts)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SCBRAW ...... Hit any key to continue ......')
pause
[As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scbraw(A,B,C,D)
err=norm(As+Bt(:,1:rank(D))*inv(Dt(1:rank(D),1:rank(D)))*Ct(1:rank(D),:)-inv(Gms)*A*Gms)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SCB ...... Hit any key to continue ......')
pause
[As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scb(A,B,C,D)
err=norm(As+Bt(:,1:rank(D))*Ct(1:rank(D),:)-inv(Gms)*A*Gms)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DSCB ...... Hit any key to continue ......')
pause
[As,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = dscb(A,B,C,D)
eerr=norm(As+Bt(:,1:rank(D))*Ct(1:rank(D),:)-inv(Gms)*A*Gms)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function KCF ...... Hit any key to continue ......')
pause
[Ks,U,V,dims] = kcf(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function MORSEIDX ...... Hit any key to continue ......')
pause
[I1,I2,I3,I4] = morseidx(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function BLKZ ...... Hit any key to continue ......')
pause
bzero = blkz(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function INVZ ...... Hit any key to continue ......')
pause
[invzs, I1] = invz(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function INFZ ...... Hit any key to continue ......')
pause
infzs = infz(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function L_INVT ...... Hit any key to continue ......')
pause
lefts = l_invt(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function R_INVT ...... Hit any key to continue ......')
pause
rights = r_invt(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function NORMRANK ...... Hit any key to continue ......')
pause
NR = normrank(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function V_STAR ...... Hit any key to continue ......')
pause
V = v_star(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function V_MINUS ...... Hit any key to continue ......')
pause
V = v_minus(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function V_PLUS ...... Hit any key to continue ......')
pause
V = v_plus(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function S_STAR ...... Hit any key to continue ......')
pause
S = s_star(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function S_MINUS ...... Hit any key to continue ......')
pause
S = s_minus(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function S_PLUS ...... Hit any key to continue ......')
pause
S = s_plus(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function R_STAR ...... Hit any key to continue ......')
pause
R = r_star(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function N_STAR ...... Hit any key to continue ......')
pause
N = n_star(A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function S_LAMBDA ...... Hit any key to continue ......')
pause
S = s_lambda(A,B,C,D,rand)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function V_LAMBDA ...... Hit any key to continue ......')
pause
V = v_lambda(A,B,C,D,rand)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SSORDER ...... Hit any key to continue ......')
pause
ss = ssorder(S,V)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SSINTSEC ...... Hit any key to continue ......')
pause
ss = ssintsec(S,V)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SSADD ...... Hit any key to continue ......')
pause
ss = ssadd(S,V)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function EA_DS ...... Hit any key to continue ......')
pause
E=B*B'*C'*C
rank_E = rank(E)
pause
[Et,At,P,Q,n1,n2] = ea_ds(E,A)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SD_DS ...... Hit any key to continue ......')
pause
[Es,As,Bs,Cs,Ds,Ez,Psi,Psc,Psd,Psr,Gme,Gms,Gmo,iGmi,dim]= sd_ds(E,A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function INVZ_DS ...... Hit any key to continue ......')
pause
[zrs,I1] = invz_ds(E,A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function INFZ_DS ...... Hit any key to continue ......')
pause
infzs = infz_ds(E,A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function L_INVT_DS ...... Hit any key to continue ......')
pause
lefts = l_invt_ds(E,A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function R_INVT_DS ...... Hit any key to continue ......')
pause
rights = r_invt_ds(E,A,B,C,D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function MPFACT ...... Hit any key to continue ......')
pause
[Am,Bm,Cm,Dm,Av,Bv,Cv,Dv] = mpfact(A,B,C,D)
s=rand;
Iv=eye(max(size(Av)));
Gm_zrs=invz(Am,Bm,Cm,Dm)
Vs2I=[Cv*(s*Iv-Av)^(-1)*Bv+Dv]*[Cv*(-s*Iv-Av)^(-1)*Bv+Dv]'
Gs_GmV=C*(s*eye(n)-A)^(-1)*B+D-[Cm*(s*eye(n)-Am)^(-1)*Bm+Dm]*[Cv*(s*Iv-Av)^(-1)*Bv+Dv]
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function IOFACT ...... Hit any key to continue ......')
pause
A=A-(max(abs(eig(A)))+1)*eye(size(A,1));
[Ai,Bi,Ci,Di,Ao,Bo,Co,Do] = iofact(A,B,C,D)
s=rand;
Ii=eye(max(size(Ai)));
Go_zrs=invz(Ao,Bo,Co,Do)
Gi2I=[Ci*(-s*Ii-Ai)^(-1)*Bi+Di]'*[Ci*(s*Ii-Ai)^(-1)*Bi+Di]
Gs_GiGo=C*(s*eye(n)-A)^(-1)*B+D-[Ci*(s*Ii-Ai)^(-1)*Bi+Di]*[Co*(s*eye(n)-Ao)^(-1)*Bo+Do]
disp('   ')
disp('...... Hit any key to continue ......')
A=-A;
pause

clc
disp('   ')
disp('...... Testing m-Function GCFACT ...... Hit any key to continue ......')
pause
[Am,Bm,Cm,Dm,Au,Bu,Cu,Du] = gcfact(A,B,C,D)
s=rand;
Iu=eye(max(size(Au)));
Gm_zrs=invz(Am,Bm,Cm,Dm)
Gs_GmU=C*(s*eye(n)-A)^(-1)*B+D-[Cm*(s*eye(n)-Am)^(-1)*Bm+Dm]*[Cu*(s*Iu-Au)^(-1)*Bu+Du]
disp('   ')
disp('...... Hit any key to continue ......')
pause

%clc
%disp('   ')
%disp('...... Testing m-Function DMPFACT ...... Hit any key to continue ......')
%pause
%[Am,Bm,Cm,Dm,Av,Bv,Cv,Dv] = dmpfact(A,B,C,D)
%z=rand;
%Iv=eye(max(size(Av)));
%Gm_zrs=invz(Am,Bm,Cm,Dm)
%Vz2I=[Cv*(z*Iv-Av)^(-1)*Bv+Dv]*[Cv*(1/z*Iv-Av)^(-1)*Bv+Dv]'
%Gz_GmV=C*(z*eye(n)-A)^(-1)*B+D-[Cm*(z*eye(n)-Am)^(-1)*Bm+Dm]*[Cv*(z*Iv-Av)^(-1)*Bv+Dv]
%disp('   ')
%disp('...... Hit any key to continue ......')
%pause

%clc
%disp('   ')
%disp('...... Testing m-Function DIOFACT ...... Hit any key to continue ......')
%pause
%A=A/(max(abs(eig(A)))+1);
%[Ai,Bi,Ci,Di,Ao,Bo,Co,Do] = diofact(A,B,C,D)
%z=rand;
%Ii=eye(max(size(Ai)));
%Go_zrs=invz(Ao,Bo,Co,Do)
%Gi2I=[Ci*(1/z*Ii-Ai)^(-1)*Bi+Di]'*[Ci*(z*Ii-Ai)^(-1)*Bi+Di]
%Gs_GiGo=C*(z*eye(n)-A)^(-1)*B+D-[Ci*(z*Ii-Ai)^(-1)*Bi+Di]*[Co*(z*eye(n)-Ao)^(-1)*Bo+Do]
%disp('   ')
%disp('...... Hit any key to continue ......')
%pause

clc
disp('   ')
disp('...... Testing m-Function SA_SEN ...... Hit any key to continue ......')
pause
c = sa_sen(A,B)
[At,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scb(A,B,c,c*B*0)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function SA_ACT ...... Hit any key to continue ......')
pause
b = sa_act(A,C)
[At,Bt,Ct,Dt,Gms,Gmo,Gmi,dim] = scb(A,b,C,C*b*0)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function ATEA ...... Hit any key to continue ......')
pause
x=rand;
if x>0.5
    Op=1;
else
    Op=0;
end
F=atea(A,B,C,D)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... The following tests require a disturbance matrix E ......')
disp('   ')
chs=input('...... Choose 1 to enter the matrix. Otherwise, it will be done by computer. ( 1 or 0 )? ');
if chs==1
    disp('   ')
    E=input('...... Please enter matrix E now. E = ');
else
    E=rand(n,1)
end
pause

clc
disp('   ')
disp('...... Testing m-Function GM2STAR ...... Hit any key to continue ......')
pause
gms2 = gm2star(A,B,C,D,E)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H2CARE ...... Hit any key to continue ......')
pause
dd=rand(size(D));
[P,err] = h2care(A,B,C,dd)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H2STATE ...... Hit any key to continue ......')
pause
F = h2state(A,B,C,D,E)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function GM8STAR ...... Hit any key to continue ......')
pause
gms8 = gm8star(A,B,C,D,E)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H8CARE ...... Hit any key to continue ......')
pause
disp('   ')
gamma=input('...... Enter a positive value for gamma = ');
disp('   ')
[P,err] = h8care(A,B,C,rand(size(D)),E,gamma)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H8STATE ...... Hit any key to continue ......')
pause
F = h8state(A,B,C,D,E,gamma)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function ADDPS ...... Hit any key to continue ......')
pause
F = addps(A,B,C,D,E,Op)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DATEA ...... Hit any key to continue ......')
pause
F = datea(A,B,C,D)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DARE ...... Hit any key to continue ......')
pause
[P,err] = dare(A,B,C'*C,D'*D,C'*D)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DGM2STAR ...... Hit any key to continue ......')
pause
gms2 = dgm2star(A,B,C,D,E)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H2DARE ...... Hit any key to continue ......')
pause
P = h2dare(A,B,C,dd)
err=norm(A'*P*A+C'*C-(A'*P*B+C'*dd)*(dd'*dd+B'*P*B)^(-1)*(A'*P*B+C'*dd)'-P)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DH2STATE ...... Hit any key to continue ......')
pause
F = dh2state(A,B,C,D,E)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DGM8STAR ...... Hit any key to continue ......')
pause
gms8 = dgm8star(A,B,C,D,E)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function H8DARE ...... Hit any key to continue ......')
pause
disp('   ')
gamma=input('...... Enter a positive value for gamma = ');
disp('   ')
P = h8dare(A,B,C,dd,E,gamma)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DH8STATE ...... Hit any key to continue ......')
pause
F = dh8state(A,B,C,D,E,gamma)
poles=eig(A+B*F)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DADDPS ...... Hit any key to continue ......')
pause
F = daddps(A,B,C,D,E)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function DDPCM ...... Hit any key to continue ......')
disp('   ')
D1=C1*E*0;
D22=C2*E*0;
disp('...... Enter system data A, B, E, C1, D1, C2, D2, D22 and issue RETURN Command ......')
pause
K = ddpcm(A,B,E,C1,D1,C2,D2,D22)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Testing m-Function ROSYS4DDP ...... Hit any key to continue ......')
pause
[Ar,Br,Er,C1r,D1r,C2r,D2r] = rosys4ddp(A,B,E,C1,D1,C2,D2,D22)
disp('   ')
disp('...... Hit any key to continue ......')
pause

clc
disp('   ')
disp('...... Congratulation! This round of test is successful! ......')