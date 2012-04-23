!********************************************
! Command Line Input Generator for ADF 
! author: Holger Kruse (holger.kruse@uni-muenster.de)
! idea and concept (and some code..) based on 'cefine' by Stefan Grimme 
! Last change: 1.4.10
! To Do:
! UKS setup 
! ciadfrc: more options?
! easier symmetry input 
!********************************************
      Program comand_line_define_adf
      implicit none
      
      integer maxarg,io,i,nn
      parameter (maxarg=40)
      character*80 arg(maxarg)    
      character*80 outfile,run1,run2,run3
      character*80 func ,fcore          
      character*80 bas,home,optim           
      character*80 sym,jobfile,pwd            
      character*80 atmp,cpid           
      character*80 infile,fragfile    
      character*80 solvent
      character*10 version
      character*250 string,sumFrag,string2,fraglist(10)
      logical KEEP,HYBRID,metahybrid,FC,GGA,SETSYM,META,ZORA,set_TC
      logical da,FON,TS,COSMO,OPTI,ECHO,RANGST,JOB,AFREQ,RESTART,TE
      logical strange_elem,COORD,VDW,PR,OPT,NOVDW,FXYZ,FTMOL,SMOOTH
      logical LWALL,FRAG,QUICK,SFREQ,QUILD,LDA,POST_META,POST_HF,SMEAR
      logical DFTB,ANC,TAPE,CNEW
      integer charge,l,ntypes,scfcv,scfcv2,nheavy
      integer nat,gcart,err,ierr,mode,pid,wall,nfrag
      real*8  desythr,xx(5),thize,dum,grid,grid2,trustrad
      real*8 xyz(3,10000)
      integer iat(10000),Slow,Shigh,giter,siter
      version='1.2'
      io=1
      
      LDA=.false.
      GGA=.true.   ! assuming always GGA
      META=.false.
      HYBRID=.false.
      METAHYBRID=.false.
      ECHO =.false.

ccccccccccccccccccccccccccccccccccc
cccccccccc   defaults ccccccccccccc
ccccccccccccccccccccccccccccccccccc
      pid=getpid()
      write(atmp,'(i7)') pid
      cpid=trim(adjustl(atmp))
      jobfile='job.'//cpid
      outfile=trim(jobfile)  
      infile='XYZ.in'
      VDW  =.false.
      bas    ='TZ2P'
      charge =0
      func ='BP86'
      grid = 7.d0
      gcart=4
      scfcv=6
      scfcv2=6
      COSMO=.false. 
      solvent='water'
      fcore='none'
      setsym=.false.
      trustrad=0.3
      optim='Delocal'
      FTMOL=.false.
      JOB=.false.
      AFREQ=.false.
      TS=.false.
      ZORA=.false.
      RESTART=.true.
      wall=120
      LWALL=.false.
      OPT=.false.
      FRAG=.false.
      fragfile='frag'
      QUICK=.false.
      SFREQ=.false.
      giter=250
      siter=125
      POST_META=.false.
      POST_HF=.false.
      SMOOTH=.false.
      set_TC=.true.
      TE=.false.
      SMEAR=.true.
      DFTB=.false.
      QUILD=.false.
      ANC=.false. 
      TAPE=.false.
      CNEW=.true.
ccccccccccccccccccccccccc
       write(*,*)'  CIADF  V ', trim(version)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc    read and evaluate /$HOME/.ciadfrc  file c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      call get_environment_variable('PWD',pwd)
      call get_environment_variable('HOME',atmp)
      home=trim(atmp)//'/.ciadfrc'
!      write(*,*) home
      inquire(file=home,exist=da)
!      inquire(file='~/.ciadfrc',exist=da) ! doesnt work with gfortran
      if(da)then
         open(unit=20,file=home)
!         open(unit=20,file='~/.ciadfrc')   ! doesnt work with gfortran
 842     read(20,'(a)',end=942)atmp
         if(index(atmp,'func').ne.0)then         
            call backstring(atmp,func,4)
         endif
         if(index(atmp,'basis').ne.0)then         
            call backstring(atmp,bas,5)
         endif
         if(index(atmp,'fc ').ne.0)then         
            call backstring(atmp,fcore,2)
         endif
         if(index(atmp,'notc').ne.0)then         
            set_TC=.false.
         endif
         if(index(atmp,'smear').ne.0)then         
            smear=.true.
         endif
         if(index(atmp,'tenergy').ne.0)then         
            TE=.true.
         endif
         if(index(atmp,'vdw').ne.0) then
            if(index(atmp,'on').ne.0)VDW=.true.   
         endif
         if(index(atmp,'trustrad').ne.0)then         
            call readl(atmp,xx,nn)
            trustrad=xx(nn)
         endif
         if(index(atmp,'grid').ne.0)then         
            call readl(atmp,xx,nn)
            grid=xx(nn)
         endif
         if(index(atmp,'maxcycle').ne.0)then         
            call readl(atmp,xx,nn)
            giter=xx(nn)
         endif
         if(index(atmp,'scfiter').ne.0)then         
            call readl(atmp,xx,nn)
            siter=xx(nn)
         endif
         if(index(atmp,'format').ne.0)then         
           if(index(atmp,'xzy').ne.0)then         
            FXYZ=.true.
           elseif(index(atmp,'tmol').ne.0)then         
            FTMOL=.true.
           endif
         endif
         if(index(atmp,'cosmo').ne.0)then         
            call backstring(atmp,solvent,5)
         endif
         if(index(atmp,'restart').ne.0)then         
            RESTART=.true.
           if(index(atmp,'off').ne.0)then
           RESTART=.false.
           endif         
         endif
         if(index(atmp,'ifile').ne.0)then         
            call backstring(atmp,infile,5)
         endif
         if(index(atmp,'ofile').ne.0)then         
            call backstring(atmp,outfile,5)
         endif
         if(index(atmp,'ADFPREP').ne.0) JOB=.true.         
         goto 842
 942     close(20)
      endif

cccccccccccccccccccccccccccccccccccccccccccc
cccccccccc   get command line arguments   cc
cccccccccccccccccccccccccccccccccccccccccccc
      arg=''
      do i=1,maxarg
         call getarg(i,arg(i))
      enddo
ccc first print help
      if(arg(1).eq.'-h' .or. arg(1).eq.'?' .or.
     .   arg(1).eq.'-help')then
      pr=.true. 
      if(pr) then
        write(*,*) '                                       '
      write(*,'(10x,''*** Command line Input generator for ADF***'')') 
        write(*,'(10x,''          H. Kruse Feb. 2010      '')') 
      write(*,'(10x,''            | V '',a,'' |       '')')trim(version)
        write(*,'(10x,''          developer version     '')') 
      endif
         write(*,*)'OPTIONS:'
         write(*,*)'   -v (print defaults, after reading .ciadfrc)'
         write(*,*)'   -h / -help / -h (this help output)'
         write(*,*)' more information in the manual'
         write(*,*)'   '
         write(*,*)'* input/ouput:'
         write(*,*)'   -n (name eg. H2_opt, submit this)'
         write(*,*)'   -f   <string>  input file (def: XYZ.in)'
         write(*,*)'   -tm (TMOL input format)'
         write(*,*)'   -xyz (XYZ (xmol) input format (def))'
         write(*,*)'  [ -job (sets JOB option=true)                ]'
         write(*,*)'  [ -o    <string> actual JOB file)            ]' 
         write(*,*)'  [ -wall <int>  (walltime in h for "jobfile") ]'
         write(*,*)' needs: <input> file in XYZ or TMOL format'
         write(*,*)' '
         write(*,*)'* calculation setup:'
         write(*,*)'   -func   <string>'
         write(*,*)'   -bas    <string>'
         write(*,*)'   -fc     <string> (frozen core, def: none)'
         write(*,*)'   -grid   <integer>'
         write(*,*)'   -vdw (DFT-D) // -novdw'
         write(*,*)'   -noopt / -nopt (do single point)'
         write(*,*)'   -opt  (do optimization)'
         write(*,*)'   -chrg   <integer>'
         write(*,*)'   -scfcv  <integer> (scf convergence, def: 7)'
         write(*,*)'   -siter  <integer> (max scf iterations)'
         write(*,*)'   -grad   <integer> (gradient convergence,def: 4)'
         write(*,*)'   -c      <integer> (max cycle in GeoOpt)'
         write(*,*)'   -cart (cartesian coordinates, def: delocal)'
         write(*,*)'   -intal (internal coordinates)'
         write(*,*)'   -sym    <string> (symmetry, eg: "C\(4h\)" )'
         write(*,*)'   -cosmo  <string> (COSMO with solvent=<string>)' 
         write(*,*)'   -afreq (analytical freq)'
         write(*,*)'   -sfreq  <int> <int> (scan freq: low high)'
         write(*,*)'   -pmeta / -phf  (POST-SCF options)'
         write(*,*)'   -zora (turns scalar zora on)'
         write(*,*)'   -smooth  (smoothing of the gradient)' 
         write(*,*)'   -ts <integer> (TS search, <int>=modes)'
         write(*,*)'   -tape (saves TAPE21 file)'
         write(*,*)'   '
         write(*,*)'* special calls:'
         write(*,*)'   -frag (handle self defined fragment files)'
         write(*,*)'   -quick (pbe-d/dz grid=3,cart,Geocycle=30)'
         write(*,*)'   -dftb (dftb input, quick pre-opt,-noopt for SP)'
         write(*,*)'   -quild (quild input with deloc coordinates)'
         write(*,*)'   -ancopt (prepares ancopt input)'
         write(*,*)'   '
         write(*,*)'   '
         stop
      endif
c now process the arguments. (A 'space' after the option is a good idea)
      do i=1,maxarg
         if(arg(i).ne.'')then
            if(index(arg(i),'-vdw').ne.0)   VDW=.true. ! DFT-D on
            if(index(arg(i),'-tape').ne.0)   TAPE=.true. ! save tape21 file
            if(index(arg(i),'-ancopt ').ne.0) then  ! ancopt external optimizer
            ANC=.true.  
            JOB=.true.
            giter=1
            optim='cartesian'
            gcart=8  ! prevents ADF from convergence
            opt=.true.
            endif
            if(index(arg(i),'-v ').ne.0)   ECHO=.true. ! print some defaults
            if(index(arg(i),'-quick ').ne.0)   QUICK=.true. ! quick&dirty opt
            if(index(arg(i),'-dftb').ne.0)   DFTB=.true. ! DFTB
            if(index(arg(i),'-quild ').ne.0)   QUILD=.true. ! quild
            if(index(arg(i),'-n ').ne.0) then       ! name of the jobfile
            jobfile=arg(i+1) 
            endif
            if(index(arg(i),'-afreq').ne.0)   AFREQ=.true. ! analytical freq
            if(index(arg(i),'-noopt').ne.0)   OPT=.false. ! SP
            if(index(arg(i),'-nopt').ne.0)   OPT=.false. ! SP
            if(index(arg(i),'-opt').ne.0)   OPT=.true. ! GEO OPT
            if(index(arg(i),'-pmeta ').ne.0)   POST_META=.true. ! meta-ggas in post-scf manner
            if(index(arg(i),'-phf ').ne.0)   POST_HF=.true. ! hybrids in pst scf manner
            if(index(arg(i),'-novdw').ne.0) VDW=.false.  ! DFT-D off
            if(index(arg(i),'-zora').ne.0) ZORA=.true.  ! scalar zora 
            if(index(arg(i),'-smooth').ne.0) SMOOTH=.true.  ! smoothing of gradient 
            if(index(arg(i),'-njob').ne.0) JOB=.false.  ! smoothing of gradient 
            if(index(arg(i),'-frag').ne.0) then  ! provide own fragments
              FRAG=.true.  ! fragments 
              fragfile=arg(i+1)
            endif
            if(index(arg(i),'-xyz ').ne.0) then ! xzy files as input
             FXYZ=.true.
             FTMOL=.false.
            endif
            if(index(arg(i),'-norestart ').ne.0)then  !infile
               restart=.false.           
            endif
            if(index(arg(i),'-sfreq ').ne.0)then  !scanfreq 
               call readl(arg(i+1),xx,nn)
               Slow=idint(xx(1))
               call readl(arg(i+2),xx,nn)
               Shigh=idint(xx(1))
               SFREQ=.true.           
            endif
            if(index(arg(i),'-job ').ne.0)then  !job
               JOB=.true.           
            endif
            if(index(arg(i),'-tm ').ne.0) then ! turbomole format as input
            FTMOL=.true.
            FXYZ=.false.
            endif
            if(index(arg(i),'-f ').ne.0)then  !infile
               infile=arg(i+1)           
            endif
            if(index(arg(i),'-wall ').ne.0)then  ! copy back important files bevor scratch is erases
               call readl(arg(i+1),xx,nn)
               wall=idint(xx(1))
               lwall=.true.        
            endif
            if(index(arg(i),'-o ').ne.0)then ! OUTPUT FILE
               outfile=arg(i+1)           
            endif
            if(index(arg(i),'-chrg').ne.0)then !Charge
               call readl(arg(i+1),xx,nn)
               charge=idint(xx(1))
            endif
            if(index(arg(i),'-cosmo').ne.0)then ! COSMO
               COSMO=.true.
               solvent=arg(i+1)
            endif
            if(index(arg(i),'-ts ').ne.0)then ! TS search
               call readl(arg(i+1),xx,nn)
               TS=.true.
               mode=idint(xx(1))
            endif
            if(index(arg(i),'-scfcv').ne.0)then ! SCFCONV
               call readl(arg(i+1),xx,nn)
               scfcv=idint(xx(1))
            endif
            if(index(arg(i),'-grad').ne.0)then ! GRADCONV
               call readl(arg(i+1),xx,nn)
               gcart=idint(xx(1))
            endif
            if(index(arg(i),'-siter ').ne.0)then ! GRADCONV
               call readl(arg(i+1),xx,nn)
               siter=idint(xx(1))
            endif
            if(index(arg(i),'-c ').ne.0)then ! geometry cycles
               call readl(arg(i+1),xx,nn)
               giter=idint(xx(1))
            endif
            if(index(arg(i),'-sym').ne.0)then  ! SET SYMMETRY
            SETSYM=.true.
               sym=arg(i+1)
            endif
            if(index(arg(i),'-grid').ne.0) then ! NUM INT Grid
             call readl(arg(i+1),xx,nn)
             grid =xx(1)
c             write(*,*) grid, xx(1)
            endif
            if(index(arg(i),'-bas').ne.0) then
             bas  =arg(i+1)  ! basis set
            endif
            if(index(arg(i),'-cart').ne.0) then ! cartesian optimizer
             optim  ='Cartesian'
            endif
            if(index(arg(i),'-intal').ne.0) then ! internal/ Zmatrix optimizer
             optim  ='Internal'
             write(*,*) 'Adjust ATOMS BLOCK by hand!!'  ! user should provide good z-matrix
            endif
            if(index(arg(i),'-fc').ne.0) then ! frozen core
             fcore  =arg(i+1)  
            endif
            if(index(arg(i),'-func').ne.0) then           ! functional
               func=arg(i+1)           
            endif
            if(index(arg(i),'-angst').ne.0)RANGST=.true. ! change cartesian coords unit to angstrom
c keep outputs for debuging purposes  (no purpose yet)
            if(index(arg(i),'-keep').ne.0)KEEP=.true. 
         endif
      enddo
ccccccccccccccccccccccccc 
      If(.not.JOB) outfile=trim(jobfile)


ccccccccccccccccccccccccccc
cccc read coordinates     c
ccccccccccccccccccccccccccc
c  more formats?
!       if(.not.Frag)
       if (FXYZ.and..not.frag)call xyzrd(xyz,iat,nat,infile)
       if (FTMOL)call tmolrd(xyz,iat,nat,infile)
ccccccccccccccccccccccccccccccccccc 
ccc prepare coord file for ANCOPT c
ccccccccccccccccccccccccccccccccccc
        if(ANC) then
        open(unit=99,file='coord')
        call wtm(xyz,iat,nat,99)
        close(99)
        write(*,*)' set "$symmetry xx" in "control" file, eg xx=c4'
c        call system('touch control')
        inquire(file='control',exist=da)
        if(.not.da) call system('echo ''$symmetry xx'' > control ')
        endif

cccccccccccccccccccccccccccccccccccccc
ccccccc catch DFTB call             cc
cccccccccccccccccccccccccccccccccccccc
      if(DFTB) then
          open(unit=55,file=outfile)
          write(55,'(''$ADFBIN/dftb <<eor '')')
          write(55,'(''ATOMS'')')
          call wxyz(xyz,iat,nat,55)
          write(55,'(''END'')')
          write(55,'(''Geometry'')')
          write(55,'(''RunType GO'')')
          write(55,'(''Optim Cartesian'')')
          write(55,'(''Iterations 10'')')
          write(55,'(''Step Trustradius=0.3'')')
          write(55,'(''END'')')
          write(55,'(''eor'')')
          close(55)
          write(*,'(5x,''output : '',a)') trim(outfile)
          stop ': DFTB quick input done'
       endif
ccccccccccccccccccccccccccccccc
ccccc   ECHO DEFAULTS        cc
ccccccccccccccccccccccccccccccc
      if(ECHO) then
      write(*,*) ' '
      write(*,*) 'some default settings (after reading .ciadfrc) '
      write(*,*)'functional  ',trim(func), ', basis set ', trim(bas) 
      write(*,'('' num grid '',f4.1,'' grad conv'',i2)')grid,gcart
      write(*,*)'frozen core ',trim(fcore), ', vdw ', VDW
      write(*,*)'COSMO ',COSMO,', solvent ', trim(solvent)
      write(*,*)'inputfile ', trim(infile),', outputfile ',trim(outfile)
      write(*,'('' scfconv '',i2,'', optim '',a)')scfcv,trim(optim)
      stop 
      endif

ccccccccccccccccccccccccccccccccccccc     
cccccccccc   QUICK Option           c
ccccccccccccccccccccccccccccccccccccc
       if(QUICK) then
       grid=4
       func='pbe'
       fcore='small'
       bas='DZ'
       scfcv=6
       gcart=3
       vdw=.true.
       optim='Cartesian'
       endif
     

      grid2=grid-1        ! whats best??
ccccccccccccccccccccccccccccccccccccc
ccccccc Print what you are doing   cc
ccccccccccccccccccccccccccccccccccccc
          if(vdw) then
           write(*,'(5x,a,''/'',a,"-D")') trim(func),trim(bas)
          else
           write(*,'(5x,a,''/'',a)') trim(func),trim(bas)
          endif
       write(*,'(5x,''INT GRID='',f5.1,'' FCORE='',a)') grid,trim(fcore)
        if(setsym) write(*,'(5x,''explixit symmetry='',a)') trim(sym)

cccccccccccccccccccccccccccccccccccc
cccc check the functional type  cccc 
cccccccccccccccccccccccccccccccccccc
c hybrids
      if( index(func,'b3lyp').ne.0 
     |.or.index(func,'bhandhlyp').ne.0 
     |.or.index(func,'b1lyp' ).ne.0 
     |.or.index(func,'b1pw91' ).ne.0 
     |.or.index(func,'o3lyp' ).ne.0 
     |.or.index(func,'x3lyp' ).ne.0 
     |.or.index(func,'kmlyp' ).ne.0 
     |.or.index(func,'mpw1pw' ).ne.0 
     |.or.index(func,'mpw1k' ).ne.0 
     |.or.index(func,'opbe0' ).ne.0 
     |.or.index(func,'pbe0'  ).ne.0 
     | ) then
         HYBRID=.true.
         GGA=.false.
      endif
c check meta-hybrid
      if( index(func,'m06-2x').ne.0
     |.or.index(func,'m06 ').ne.0
     |.or.index(func,'m06-hf' ).ne.0
     |.or.index(func,'tpssh'  ).ne.0
     | ) then
         METAHYBRID=.true.
         GGA=.false.
      endif    
c check meta-GGA
      if( index(func,'tpss ').ne.0 
     |.or.index(func,'m06-l').ne.0 
     |.or.index(func,'m06l').ne.0 
     |.or.index(func,'ssb-d' ).ne.0 
     | ) then
       if (grid.lt.6) then
       write(*,*) 'Enforcing integration grid 6 (meta-GGA)!'   !numerically unstable otherwise
       grid=6
       endif
         META=.true.
         GGA=.false.
      endif      
      if( index(func,'vwn').ne.0 
     |.or.index(func,'pw92').ne.0 
     |.or.index(func,'lda' ).ne.0 
     |.or.index(func,'xalpha'  ).ne.0 
     | ) then
      GGA=.false.
      LDA=.true.
      endif
      if(HYBRID) write(*,*)'hybrid functional detected'
      if(META) write(*,*)'meta GGA functional detected'
      if(METAHYBRID) write(*,*)'meta hybrid functional detected'
      if(LDA) write(*,*)'LDA functional detected'
cccccccccccccccccccccccccc
c check special cases:   c
cccccccccccccccccccccccccc
      if (index(func,'ssb-d').ne.0) then
      VDW=.false.
      endif
cccccccccccccccccccccccccccc
c check for user errors:   c
cccccccccccccccccccccccccccc
      !metaGGA/ F_exchange and frozen core are no good 
       if(fcore.ne.'none'.and.(POST_META.or.META))
     . stop ': frozen core and metaGGAs are no good'
       if(fcore.ne.'none'.and.(POST_HF.or.HYBRID))
     . stop ': frozen core and Hybrids are no good'

      !scanfreq
      if(Slow.gt.Shigh) stop   !check for input error
     . ': interchange scan frequencies. low<->high '
ccccccccccccccccccccccccccccccccccccccccccccc

! handle fragment files, fragment file name given by the user
! expecting: <Fragment name>.t21 and <fragment name>.xyz
! <Fragment name>
! <Fragment name2>
      IF(FRAG) then
      inquire(file=fragfile,exist=da)
      if(da) then
! read fragfile
        sumFrag=''
        nFrag=0
        io=2
        open(unit=io,file=fragfile)
 142    read(io,'(a)',end=242)atmp   
        nfrag=nfrag+1
        fraglist(nfrag)=trim(atmp)
        goto 142
 242    close(io)
      else
        write(*,*)'No fragment file found!'
        write(*,*)' handle input yourself '
        FRAG=.false.
      endif
      write(*,'(5x,''Found '',i2,'' fragments !'')') nfrag
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc   print output file                 cccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      inquire(file=outfile,exist=da)
      if(da) then
      call system('cp '//outfile//' input.bak')
       write(*,*) '  !!!  Warning input file already there  !!!'
       write(*,*) '  old input now: input.bak'
      endif
      if (JOB) outfile=trim(jobfile)//'_'//trim(outfile)
      if(ANC) outfile='adf.in'
      if(io.ne.6)open(unit=io,file=outfile)
      if (set_TC) pwd='$TC_SUBMISSION_DIR'
! START PARSING 
      if (.not.JOB) then

        if (.not.quild) then
        write(io,
     .  '(''$ADFBIN/adf <<eor>'',a,''/''a,''.out'')')
     .  trim(pwd) ,trim(jobfile)
        else
         write(io,
     .  '(''$ADFBIN/quild <<eor>'',a,''/''a,''.out'')')
     .  trim(pwd) ,trim(jobfile)
      endif
      endif
      write(io,'(''TITEL '',a)') trim(jobfile) !titel
      write(io,'('' '')') !free line



! XC BLOCK
      write(io,'(''XC '')')
      if (LDA) write(io,'(''LDA  '',a)') trim(func)
      if (GGA) write(io,'(''GGA  '',a)') trim(func)
      if (META) write(io,'(''MetaGGA  '',a)') trim(func)
      if (HYBRID) write(io,'(''Hybrid  '',a)') trim(func)
      if (METAHYBRID) write(io,'(''metaHybrid  '',a)') trim(func)
      if (VDW) write(io,'(''DISPERSION'')')
      if(POST_META) write(io,'(''METAGGA '')')
      if(POST_HF) write(io,'(''HARTREEFOCK '')')
      write(io,'(''END '')')
      write(io,'('' '')')

! Relativistic
      if (ZORA) then
      write(io,'('' '')')
      write(io,'(''RELATIVISTIC Scalar ZORA'')')
      write(io,'('' '')')
      endif

! PRINT BLOCK
      write(io,'(''EPRINT '')')
      if (TAPE) write(io,'(''SFO eig ovl '')')
      if (.not.TAPE) write(io,'(''SFO noeig noovl '')')
      write(io,'(''END '')')
      write(io,'('' '')')

!BASIS BLOCK
      write(io,'(''BASIS '')')
      write(io,'(''TYPE '',a)') trim(bas)
      write(io,'(''CORE '',a)') trim(fcore)
      write(io,'(''END '')')
      write(io,'('' '')')

! NUMERICAL INTEGRATION
      write(io,'(''INTEGRATION '',2f5.1)') grid , grid2
      write(io,'('' '')')

      if(.not.quild) then
! GEOMETRY BLOCK
       write(io,'(''GEOMETRY '')')
       if(OPT) then
         write(io,'(''Optim '',a)') trim(optim)
       else
       write(io,'(''SP'')')
       endif
       if(.not.quick) write(io,'(''Iterations '',i4)') giter
       if(QUICK) write(io,'(''Iterations 30 '')')
       write(io,'(''converge 1.0e-'',i1)')gcart
       write(io,'(''Step TrustRadius='',f4.2)') trustrad
       if(TS) write(io,'(''TransitionState mode='',i1)') mode
       if(SMOOTH) write(io,'(''smooth freezecells'')') 
       write(io,'(''END '')')
       write(io,'('' '')')
      else
       write(io,'(''GEOMETRY '')')
       write(io,'(''END '')')
       write(io,'('' '')')
      endif
!ANCOPT Options
       if(ANC) then
       write(io,'(''STOPAFTER GGRADS '')')
       write(io,'('' '')')
       endif

! QUILD BLOCK
       if(quild) then
       write(io,'(''QUILD '')')
       write(io,'(''mxgeo 300 '')')
       write(io,'(''cvg_grd=1.0e-4 '')')
       write(io,'(''icreate=7 '')')  ! whats best? 10??
       write(io,'(''IDELOCAL 1'')')
      write(io,'(''RTTRUST 0.5 '')')
       write(io,'(''END'')')
      write(io,'('' '')')
       endif

! SCF BLOCK
      write(io,'(''SCF '')')
      write(io,'(''diis ok=0.01 '')')
      write(io,'(''iterations '',i5)') siter
      write(io,'(''converge 1.0e-'',i1,''  1.0e-'',i1)') scfcv,scfcv2
      write(io,'(''END '')')
      write(io,'('' '')')

! Analytical Freq 
      if (AFREQ) then     
      write(io,'(''AnalyticalFreq '')')
      write(io,'(''End '')')
      endif
! SCANFREQ
      if(SFREQ) then
      write(io,'(''ScanFreq '',i6,i6)') Slow,Shigh
      write(io,'(''End '')')
      endif

! SYMMETRY
      if(setsym) then
      write(io,'(''SYMMETRY '',a)') trim(sym)
c      write(io,'(''End '')')
      write(io,'('' '')')
      endif
! CHARGE
      write(io,'(''CHARGE '',i2)') charge
      write(io,'('' '')')

!SOLVATION BLOCK
      if(COSMO.and.CNEW) then
      write(io,'(''INLINE $HOME/cosmo_new '')')
      endif

! SOME ADDTIONAL HELP SETTINGS
      if(.not.smear) write(io,'(''OCCUPATIONS Smear=0.0 
     . ! Avoid unphysical el-searing'')')
!      write(io,'(''smear=0 '')')
!      write(io,'(''END '')')
      write(io,'('' '')')



      write(io,'(''PRINT LOGFILE '')')
!      write(io,'(''SAVE TAPE21 TAPE13 '')')
      if (TE) write(io,'(''TOTALENERGY'')')
      write(io,'('' '')')

! ATOMS BLOCK
      write(io,'(''ATOMS '')')
      if (FRAG) then
       do i=1,nFrag
        string=trim(fraglist(i))//'.xyz'
        call xyzrd(xyz,iat,nat,string)
        call wfrag(xyz,iat,nat,io,fraglist(i))
       enddo
      else
        call wxyz(xyz,iat,nat,io)
      endif
      write(io,'(''END '')')
      write(io,'('' '')')

! FRAGMENTS BLOCK
      if(FRAG) then
      call get_environment_variable('PWD',atmp)
      write(io,'(''FRAGMENTS '')')
       do i=1,nFrag
        write(io,'(6a)')
     . trim(fraglist(i)),' ',trim(atmp),'/',trim(fraglist(i)),'.t21'
       enddo
      write(io,'(''END '')')
      write(io,'('' '')')
      endif

! END THE INPUT AND THE FILE

      write(io,'(''! RESTART ''a,
     .''/'',a,''.t21'')')trim(pwd), trim(jobfile)
!      write(io,'(''! RESTART '',a,''.t13'')') trim(jobfile)
      write(io,'(''! RESTART ''a,
     .''/'',a,''.t13'')')trim(pwd), trim(jobfile)
      write(io,'(''END INPUT '')')
      if (.not.job) then
      write(io,'(''eor '')')
      if(TAPE) write(io,'(''cp TAPE21 '',a,''/'',a,''.t21'')')
     .trim(pwd),trim(jobfile)
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccc parsing done! Not additional or alternative files  cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      If (.not.JOB) write(*,'(5x,''output : '',a)') trim(outfile)
      if(io.ne.6)close(io)
! prepare alternative jobscript, ready for submission
      if(JOB) then
      io=2
      open(unit=io,file=trim(jobfile))
!      jobfile=trim(jobfile)//'_'//trim(outfile)
      write(*,'(5x,''output : '',a)') trim(jobfile)
!      call get_environment_variable('PWD',pwd)
!      write(io,'(''PP='',a)') trim(pwd)
      write(io,'(''#!/bin/bash'')') 
      write(io,'(''cp '',a,''/'',a,'' .'')') trim(pwd),trim(outfile)
      if(RESTART) then
      write(io,'(''cp '',a,''/'',a,''.t21 .'')') trim(pwd),trim(jobfile)
!      write(io,'(''cp '',a,''/'',a,''.t13 .'')') trim(pwd),trim(jobfile)
      endif

! if your job will except the walltime and you need some files from scratch 
      if(LWALL) then
!      write(io,'(''$ADFBIN/adf < '',a,'' > '',a,''/'',a,''.out &'')')
!     . trim(outfile), trim(pwd), trim(jobfile)
      write(atmp,*)wall*60*59
!      write(io,*)'max=',wall*60*59
      write(io,'(''max='',a)')trim(adjustl(atmp)) 
      write(io,*)'nohup ./.runjob > out.tmp &'
      write(io,*)'for i in `seq 1 $max`; do'
      write(io,*)'sleep 1'
      write(io,*)'if [ -e ready ]; then'
      write(io,*)'echo "done"'
      write(io,*)'exit'
      write(io,*)'fi'
      write(io,*)'done'
      write(io,*)'rm ready'
!       write(io,*)'sleep ', wall*60*59 ! wall in h, minus 1min, in sec  
      elseif(ANC) then
      write(io,'(''cp '',a,''/'',a,'' .'')') trim(pwd),'coord'
      write(io,'(''cp '',a,''/'',a,'' .'')') trim(pwd),'control'
      write(io,'(''ancopt -adf -gcart 4 >'',a,''/'',a,''.out'')')
     . trim(pwd), trim(jobfile)
      write(io,'(''cp ADF.XYZ '',a,''/.'')') trim(pwd)
      write(io,'(''cp job.last '',a,''/.'')') trim(pwd)
      else
      write(io,'(''$ADFBIN/adf < '',a,'' > '',a,''/'',a,''.out'')')
     . trim(outfile), trim(pwd), trim(jobfile)
      endif

      if(TAPE) write(io,'(''cp TAPE21 '',a,''/'',a,''.t21'')')
     . trim(pwd),trim(jobfile)
!      write(io,'(''cp TAPE13 '',a,''/'',a,''.t13'')')
!     . trim(pwd),trim(jobfile)
      write(io,'(''cd '',a)') trim(pwd)
      write(io,'(''cp '',a,''.out adf.out'')') trim(jobfile)
      endif
      close(io)   
        
!      open(unit=io,file='.run')
!      write(io,*)'#!/bin/bash'
!      write(io,*)'max=100'
!      write(io,*)'./.runjob &'
!      write(io,*)'for i in `seq 1 $max`; do'
!      write(io,*)'sleep 1'
!      write(io,*)'if [ -e ready ]; then'
!      write(op,*)'echo "done"'
!      write(op,*)'exit'
!      write(op,*)'fi'
!      write(op,*)'done'
!      write(op,*)'rm ready'
!      close(io)
!      call system('chmod +x .run')
      if(JOB.and.LWALL) then
      open(unit=io,file='.runjob')
      write(io,'(''$ADFBIN/adf < '',a,'' > '',a,''/'',a,''.out'')')
     . trim(outfile), trim(pwd), trim(jobfile)
      write(io,*)'touch ready'
      close(io)
      call system('chmod +x .runjob')
      endif
      write(*,'(5x,''done!'')') 
      end         

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC           S. Grimme's helper subroutines                   CCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c returns the number of heavy (n>10) atoms in a coord file, and the number of diff atom types
c sets cu_pd true for some special elements that need more input
      subroutine atoms(n,nat,nt,cu_pd)
      implicit none
      integer n,i,j,nn,nat
      logical cu_pd
      character*80 a80
      real*8 xx(10)
      integer na(110),nt,na2(110)

      cu_pd=.false.
      na = 0
      na2=0
      n=0
      nat=0
      j=0
      open(unit=1,file='coord')
      read(1,'(a)',end=100) a80
 10   read(1,'(a)',end=100) a80
      call readl(a80,xx,nn)
      if(index(a80,'$').ne.0)goto 100
      if(nn.eq.3)then
         nat=nat+1
         j=j+1
         call elem(a80,i)
         na2(i)=na2(i)+1
         if(i.gt.10) na(i)=na(i)+1
cts check for Cu/Pd problem
         if(i.eq.29.or.i.eq.46) cu_pd=.true.
      endif
      goto 10
100   close(1)

      n=0
      do i=1,110
         if(na(i).gt.0)n=n+1
         if(na2(i).gt.0) nt=nt+1 ! hok   nt= number of different atom types
      enddo

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ELEM(KEY1, NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) KEY1
      CHARACTER*2 ELEMNT(107),E

      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
     
      nat=0
      e='  '
      k=1
      DO J=1,len(key1)
         if (k.gt.2)exit
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            call lower(key1(j:j))
            N=ICHAR(key1(J:J))
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
      enddo

      DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

C     *****************************************************************         

      SUBROUTINE lower(AS)
      CHARACTER*1 AS
      AS=CHAR(ICHAR(AS)-ICHAR('A')+ICHAR('a'))
      END

C     *****************************************************************         

      SUBROUTINE backstring(A1,A2,lena2)
      CHARACTER*(*) A1
      CHARACTER*(*) A2
      integer n,lena2
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo
      DO J=1,lena2   
         a2(j:j)=' '
      enddo
      a1=a2
      a2='                                                            '
      n=0
      DO J=1,len(a1)
         if(a1(j:j).ne.' ')then
            n=n+1
            a2(n:n)=a1(j:j)
         endif
      enddo

      END


C     *****************************************************************         
                                                                                
      SUBROUTINE READL(A1,X,N)                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      CHARACTER*(*) A1                                                      
      DIMENSION X(*)                                                            
      I=0                                                                       
      IS=1                                                                      
  10  I=I+1                                                                     
      X(I)=READAA(A1,IS,IB,IE)                                               
      IF(IB.GT.0 .AND. IE.GT.0) THEN                                            
                                IS=IE                                           
                                GOTO 10                                         
      ENDIF                                                                     
      N=I-1                                                                     
      RETURN                                                                    
      END                                                                       
                                                                                
                                                                                
      FUNCTION READAA(A,ISTART,IEND,IEND2)                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      REAL*8 READAA                                                             
      CHARACTER*(*) A                                                      
      NINE=ICHAR('9')                                                           
      IZERO=ICHAR('0')                                                          
      MINUS=ICHAR('-')                                                          
      IDOT=ICHAR('.')                                                           
      ND=ICHAR('D')                                                             
      NE=ICHAR('E')                                                             
      IBL=ICHAR(' ')                                                            
      IEND=0                                                                    
      IEND2=0                                                                   
      IDIG=0                                                                    
      C1=0                                                                      
      C2=0                                                                      
      ONE=1.D0                                                                  
      X = 1.D0                                                                  
      NL=LEN(A) 
      DO 10 J=ISTART,NL-1                                                       
         N=ICHAR(A(J:J))                                                          
         M=ICHAR(A(J+1:J+1)) 
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20                      
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO                            
     1 .OR. M.EQ.IDOT)) GOTO 20                                                 
   10 CONTINUE                                                                  
      READAA=0.D0                                                               
      RETURN                                                                    
   20 CONTINUE                                                                  
      IEND=J                                                                    
      DO 30 I=J,NL                                                              
         N=ICHAR(A(I:I))                                                          
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C1=C1*10+N-IZERO                                                    
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN                                     
            ONE=-1.D0                                                           
         ELSEIF(N.EQ.IDOT) THEN                                                 
            GOTO 40                                                             
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   30 CONTINUE                                                                  
   40 CONTINUE                                                                  
      IDIG=0                                                                    
      DO 50 II=I+1,NL                                                           
         N=ICHAR(A(II:II))                                                         
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN                                      
            IDIG=IDIG+1                                                         
            IF (IDIG.GT.10) GOTO 60                                             
            C2=C2*10+N-IZERO                                                    
            X = X /10                                                           
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN                                    
            X=-X                                                                
         ELSE                                                                   
            GOTO 60                                                             
         ENDIF                                                                  
   50 CONTINUE                                                                  
C                                                                               
C PUT THE PIECES TOGETHER                                                       
C                                                                               
   60 CONTINUE                                                                  
      READAA= ONE * ( C1 + C2 * X)                                              
      DO 55 J=IEND,NL                                                           
         N=ICHAR(A(J:J))                                                          
         IEND2=J                                                                
         IF(N.EQ.IBL)RETURN                                                     
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57                                           
      RETURN                                                                    
                                                                                
   57 C1=0.0D0                                                                  
      ONE=1.0D0                                                                 
      DO 31 I=J+1,NL                                                            
         N=ICHAR(A(I:I))                                                          
         IEND2=I                                                                
         IF(N.EQ.IBL)GOTO 70                                                    
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO                      
         IF(N.EQ.MINUS)ONE=-1.0D0                                               
   31 CONTINUE                                                                  
   61 CONTINUE                                                                  
   70 READAA=READAA*10**(ONE*C1)                                                
      RETURN                                                                    
      END                                                                       

C     *****************************************************************         

