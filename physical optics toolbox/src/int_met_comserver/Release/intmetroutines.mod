  u?  �   k820309    �  47D3B   6.6     ��B                                                                                                                        
       C:\My_Documents\internal_metrology\diffraction\int_met_comserver\int_met_routines.f90 INTMETROUTINES          FOCUSLENS APPLYMASK APPLYMASKROTATE APPLYMASKGENERAL APPLYMASKX APPLYMASKXROUNDED APPLYMASKPOLY                                         
                                               
                                               
            @       �   @                      
                        @                      
                        @                      
                    �                          
                  � @                       '           #DX 	   #DY 
   #WAVELENGTH    #CURVATURE    #FFTW_PLAN2D_FOR    #FFTW_PLAN2D_REV    #CZTPLAN    #AMP    #X    #Y    #LPFFTWPLANLOCK    #LPREADWRITEUNFLOCK    #DUMMY            � d                     	           
                            
                         0.0        � d                     
          
                            
                         0.0        � d                               
                            
                         0.0        � d                               
                            
                         0.0         � d                                             � d                          (                  � d                           `   0      #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS               � @                       '`            #WAVEFRONTS^FFTW_PLAN1D_X_FOR    #WAVEFRONTS^FFTW_PLAN1D_Y_FOR    #WAVEFRONTS^FFTW_PLAN1D_X_INV    #WAVEFRONTS^FFTW_PLAN1D_Y_INV    #WAVEFRONTS^V_FFTW1D_X    #WAVEFRONTS^V_FFTW1D_Y            � D                                                                                     0        � D                                                                                    0        � D                                                                                    0        � D                                                                                    0       �D                                            &                                         y                                          �D                           @                &                                         y                                        �d                            �                &           &                                         y                                        �d                            �      	  
        &                                         y
                                        �d                            �      
  
        &                                         y
                                           � d                          �                                                           0        � d                                                                                     0         � d                                             �  @                       '�           #FOCLEN    #DIAM    #XDECENTER     #YDECENTER !   #ABERR "            � D                                 
           � D                               
                            
                         0.0        � D                                
                            
                         0.0        � D                     !          
                            
                         0.0         � D                      "     �         #ZERNIKEABERRATION #              �  @                  #     '�           #OPTICS_ROUTINES^REFDIAM $   #OPTICS_ROUTINES^XO %   #OPTICS_ROUTINES^YO &   #OPTICS_ROUTINES^MODENUM '   #OPTICS_ROUTINES^COEFFVALUE (   #OPTICS_ROUTINES^ISRMS )   #OPTICS_ROUTINES^MAXMODES *           � d                    $           
                            
                         0.0        � d                    %          
                            
                         0.0        � d                    &          
                            
                         0.0        � d                     '     $             p      p $       p $                    $     $                                      0        � d                    (     $   �        
  p      p $       p $                    $     $            
                         0.0        � d                     )     �                                                     짦 x          � d                     *     �                                                          0#     @                         +                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #GETWAVEFRONTSAMPLING%SIZE ,   #GETWAVEFRONTSAMPLING%ASSOCIATED -   #GETWAVEFRONTSAMPLING%PRESENT .   #BEAM /   #NX 0   #NY 1   #DX 2   #DY 3   #XO 4   #YO 5            @                     ,     SIZE          @                     -     ASSOCIATED          @                     .     PRESENT       
   @                      /          #WAVEFRONT            @                      0                @                      1                @                     2     
           @                     3     
           @                     4     
           @                     5     
   #     @                         6                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #PROPAGATE%SIZE 7   #PROPAGATE%ABS 8   #PROPAGATE%PRESENT 9   #PROPAGATE%EPSILON :   #PROPAGATE%SIGN ;   #BEAM <   #Z =   #DXO >   #DYO ?   #APPLYCURV @            @                     7     SIZE          @                     8     ABS          @                     9     PRESENT          @                     :     EPSILON          @                     ;     SIGN      
  @                      <           #WAVEFRONT          
   @                     =     
        
  @                     >     
        
  @                     ?     
        
  @                      @       #     @                         A              	   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #CLIPRECT%ABS B   #CLIPRECT%PRESENT C   #CLIPRECT%EPSILON D   #CLIPRECT%COS E   #CLIPRECT%SIN F   #BEAM G   #LENGTH H   #WIDTH I   #XC J   #YC K   #ANGLE L   #OBSC M   #TILTANGLE N   #TILTORIENTATION O            @                     B     ABS          @                     C     PRESENT          @                     D     EPSILON          @                     E     COS          @                     F     SIN       
  @                      G           #WAVEFRONT          
   @                     H     
        
   @                     I     
        
  @                     J     
        
  @                     K     
        
  @                     L     
        
  @                      M             
  @                     N     
        
  @                     O     
  #     @                         P                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #WAVEFRONTCUMSUM%SIZE Q   #WAVEFRONTCUMSUM%ABS R   #BEAM1 S   #BEAM2 T            @                     Q     SIZE          @                     R     ABS       
  @                      S           #WAVEFRONT          
   @                      T          #WAVEFRONT    #     @                         U                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #CLEARWAVEFRONT%ASSOCIATED V   #BEAM W            @                     V     ASSOCIATED       
  @                      W           #WAVEFRONT    #     @                         X                 #CLIPPOLY%POLYGON Y   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #CLIPPOLY%ABS \   #CLIPPOLY%PRESENT ]   #CLIPPOLY%EPSILON ^   #CLIPPOLY%COS _   #CLIPPOLY%SIN `   #BEAM a   #POLY b   #XC c   #YC d   #ANGLE e   #OBSC f   #TILTANGLE g   #TILTORIENTATION h             �  @                  Y     '0            #WAVEFRONTS^NUMVERTICES Z   #WAVEFRONTS^VERTEX [           � d                     Z                                                                 0       �d                   [                
        &           &                                         y
                                            @                     \     ABS          @                     ]     PRESENT          @                     ^     EPSILON          @                     _     COS          @                     `     SIN       
  @                      a           #WAVEFRONT          
   @                      b     0      #CLIPPOLY%POLYGON Y         
  @                     c     
        
  @                     d     
        
  @                     e     
        
  @                      f             
  @                     g     
        
  @                     h     
  %     @                       i                  
   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #BEAM j         
   @                      j          #WAVEFRONT    #     @                         k              
   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #GENTOPHATBEAM%PRESENT l   #TOPHAT m   #NX n   #NY o   #DX p   #DY q   #LAMBDA r   #DBEAM s   #XC t   #YC u   #POWER v            @                     l     PRESENT       
  @                      m           #WAVEFRONT          
   @                      n             
   @                      o             
   @                     p     
        
   @                     q     
        
   @                     r     
        
   @                     s     
        
  @                     t     
        
  @                     u     
        
  @                     v     
  #     @                         w                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #CLIPCIRC%ABS x   #CLIPCIRC%PRESENT y   #CLIPCIRC%EPSILON z   #CLIPCIRC%COS {   #BEAM |   #DIAM }   #XC ~   #YC    #OBSC �   #TILTANGLE �   #TILTORIENTATION �            @                     x     ABS          @                     y     PRESENT          @                     z     EPSILON          @                     {     COS       
  @                      |           #WAVEFRONT          
   @                     }     
        
  @                     ~     
        
  @                          
        
  @                      �             
  @                     �     
        
  @                     �     
  #     @                         �                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #APPLYARBITRARYMASK%SIZE �   #APPLYARBITRARYMASK%PRESENT �   #BEAM �   #MASK �   #NEGATIVE �             @                    �     SIZE           @                    �     PRESENT       
  @                      �           #WAVEFRONT          
   @                      �          #WAVEFRONT          
  @                      �       #     @                          �                  #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #BEAM �   #FOCUSLENS_F �   #FOCUSLENS_D �   #DXOUT �   #DYOUT �         
D @@                      �           #WAVEFRONT          
  @@                     �     
        
  @@                     �     
        
  @@                     �     
        
  @@                     �     
  #     @                          �              
   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #APPLYMASK%PRESENT �   #BEAM �   #LEN_R �   #LEN_T �   #OFFSET �   #DIRECTION �   #ROTANGLE �   #X_ALIGN �   #Y_ALIGN �   #TILTANGLE �   #TILTORIENTATION �             @                    �     PRESENT       
D @@                      �           #WAVEFRONT          
   @                     �     
        
   @                     �     
        
   @                     �     
        
   @                     �            1       
 @@                     �     
        
 @@                     �     
        
 @@                     �     
        
 @@                     �     
        
 @@                     �     
  #     @                          �              
   #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #APPLYMASKROTATE%PRESENT �   #BEAM �   #LEN_R �   #LEN_T �   #OFFSET �   #ROTANGLE �   #DIRECTION �   #X_ALIGN �   #Y_ALIGN �   #TILTANGLE �   #TILTORIENTATION �             @                    �     PRESENT       
D @@                      �           #WAVEFRONT          
   @                     �     
        
   @                     �     
        
   @                     �     
        
  @@                     �     
        
   @                     �            1       
 @@                     �     
        
 @@                     �     
        
 @@                     �     
        
 @@                     �     
  #     @                          �                  #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #BEAM �   #LEN_X �   #LEN_Y �   #OFF_X �   #OFF_Y �         
D @@                      �           #WAVEFRONT          
  @@                     �     
        
  @@                     �     
        
  @@                     �     
        
  @@                     �     
  #     @                          �                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #APPLYMASKX%SQRT �   #BEAM �   #RADIUS �   #LENGTH �             @                    �     SQRT       
D @@                      �           #WAVEFRONT          
   @                     �     
        
  @@                     �     
  #     @                          �                 #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #APPLYMASKXROUNDED%SQRT �   #BEAM �   #RADIUS �   #LENGTH �   #CORNER_RAD �             @                    �     SQRT       
D @@                      �           #WAVEFRONT          
   @                     �     
        
  @@                     �     
        
   @                     �     
  #     @                          �                  #WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS    #WAVEFRONT    #BEAM �   #VERTICES �   #XC �   #YC �         
D @@                      �           #WAVEFRONT          
   @  �                   �           
          & p          & p                            
  @@                     �     
        
  @@                     �     
     �   m      fn#fn $     l   b   uapp(INTMETROUTINES    �  4   J  KINDS    �  4   J  SI_UNITS    �  4   J  CONSTANTS      4   J  WAVEFRONTS     Q  4   J  OPTICS_ROUTINES !   �  4   J  UTILITY_ROUTINES    �  4   J  ERROR_EXIT %   �  �       WAVEFRONT+WAVEFRONTS (   �  w   %   WAVEFRONT%DX+WAVEFRONTS (   T  w   %   WAVEFRONT%DY+WAVEFRONTS 0   �  w   %   WAVEFRONT%WAVELENGTH+WAVEFRONTS /   B  w   %   WAVEFRONT%CURVATURE+WAVEFRONTS 5   �  8   %   WAVEFRONT%FFTW_PLAN2D_FOR+WAVEFRONTS 5   �  8   %   WAVEFRONT%FFTW_PLAN2D_REV+WAVEFRONTS -   )  f   %   WAVEFRONT%CZTPLAN+WAVEFRONTS Z   �  �      WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS+DISCRETE_TRANSFORMS=CZTFFTWPLANS q   �  u   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%FFTW_PLAN1D_X_FOR+DISCRETE_TRANSFORMS=FFTW_PLAN1D_X_FOR q   �  u   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%FFTW_PLAN1D_Y_FOR+DISCRETE_TRANSFORMS=FFTW_PLAN1D_Y_FOR q   s  u   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%FFTW_PLAN1D_X_INV+DISCRETE_TRANSFORMS=FFTW_PLAN1D_X_INV q   �  u   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%FFTW_PLAN1D_Y_INV+DISCRETE_TRANSFORMS=FFTW_PLAN1D_Y_INV c   ]	  �   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%V_FFTW1D_X+DISCRETE_TRANSFORMS=V_FFTW1D_X c   
  �   %   WAVEFRONTS^REPLACEWAVEFRONT%CZTFFTWPLANS%V_FFTW1D_Y+DISCRETE_TRANSFORMS=V_FFTW1D_Y )   �
  �   %   WAVEFRONT%AMP+WAVEFRONTS '   Y  �   %   WAVEFRONT%X+WAVEFRONTS '   �  �   %   WAVEFRONT%Y+WAVEFRONTS 4   �  u   %   WAVEFRONT%LPFFTWPLANLOCK+WAVEFRONTS 8     u   %   WAVEFRONT%LPREADWRITEUNFLOCK+WAVEFRONTS +   �  8   %   WAVEFRONT%DUMMY+WAVEFRONTS +   �  {       SIMPLELENS+OPTICS_ROUTINES 2   >  8   %   SIMPLELENS%FOCLEN+OPTICS_ROUTINES 0   v  w   %   SIMPLELENS%DIAM+OPTICS_ROUTINES 5   �  w   %   SIMPLELENS%XDECENTER+OPTICS_ROUTINES 5   d  w   %   SIMPLELENS%YDECENTER+OPTICS_ROUTINES 1   �  O   %   SIMPLELENS%ABERR+OPTICS_ROUTINES -   *  �       ZERNIKEABERRATION+WAVEFRONTS =   )  w   %   ZERNIKEABERRATION%REFDIAM+WAVEFRONTS=REFDIAM 3   �  w   %   ZERNIKEABERRATION%XO+WAVEFRONTS=XO 3     w   %   ZERNIKEABERRATION%YO+WAVEFRONTS=YO =   �  �   %   ZERNIKEABERRATION%MODENUM+WAVEFRONTS=MODENUM C   ;  �   %   ZERNIKEABERRATION%COEFFVALUE+WAVEFRONTS=COEFFVALUE 9   �  t   %   ZERNIKEABERRATION%ISRMS+WAVEFRONTS=ISRMS ?   ^  u   %   ZERNIKEABERRATION%MAXMODES+WAVEFRONTS=MAXMODES 0   �        GETWAVEFRONTSAMPLING+WAVEFRONTS :   �  1      GETWAVEFRONTSAMPLING%SIZE+WAVEFRONTS=SIZE F     7      GETWAVEFRONTSAMPLING%ASSOCIATED+WAVEFRONTS=ASSOCIATED @   T  4      GETWAVEFRONTSAMPLING%PRESENT+WAVEFRONTS=PRESENT 5   �  C   e   GETWAVEFRONTSAMPLING%BEAM+WAVEFRONTS 3   �  0   e   GETWAVEFRONTSAMPLING%NX+WAVEFRONTS 3   �  0   e   GETWAVEFRONTSAMPLING%NY+WAVEFRONTS 3   +  0   e   GETWAVEFRONTSAMPLING%DX+WAVEFRONTS 3   [  0   e   GETWAVEFRONTSAMPLING%DY+WAVEFRONTS 3   �  0   e   GETWAVEFRONTSAMPLING%XO+WAVEFRONTS 3   �  0   e   GETWAVEFRONTSAMPLING%YO+WAVEFRONTS %   �        PROPAGATE+WAVEFRONTS /   �  1      PROPAGATE%SIZE+WAVEFRONTS=SIZE -   0  0      PROPAGATE%ABS+WAVEFRONTS=ABS 5   `  4      PROPAGATE%PRESENT+WAVEFRONTS=PRESENT 5   �  4      PROPAGATE%EPSILON+WAVEFRONTS=EPSILON /   �  1      PROPAGATE%SIGN+WAVEFRONTS=SIGN *   �  C   e   PROPAGATE%BEAM+WAVEFRONTS '   <  0   e   PROPAGATE%Z+WAVEFRONTS )   l  0   e   PROPAGATE%DXO+WAVEFRONTS )   �  0   e   PROPAGATE%DYO+WAVEFRONTS /   �  0   e   PROPAGATE%APPLYCURV+WAVEFRONTS $   �  E      CLIPRECT+WAVEFRONTS ,   A  0      CLIPRECT%ABS+WAVEFRONTS=ABS 4   q  4      CLIPRECT%PRESENT+WAVEFRONTS=PRESENT 4   �  4      CLIPRECT%EPSILON+WAVEFRONTS=EPSILON ,   �  0      CLIPRECT%COS+WAVEFRONTS=COS ,   	  0      CLIPRECT%SIN+WAVEFRONTS=SIN )   9  C   e   CLIPRECT%BEAM+WAVEFRONTS +   |  0   e   CLIPRECT%LENGTH+WAVEFRONTS *   �  0   e   CLIPRECT%WIDTH+WAVEFRONTS '   �  0   e   CLIPRECT%XC+WAVEFRONTS '     0   e   CLIPRECT%YC+WAVEFRONTS *   <  0   e   CLIPRECT%ANGLE+WAVEFRONTS )   l  0   e   CLIPRECT%OBSC+WAVEFRONTS .   �  0   e   CLIPRECT%TILTANGLE+WAVEFRONTS 4   �  0   e   CLIPRECT%TILTORIENTATION+WAVEFRONTS +   �  �       WAVEFRONTCUMSUM+WAVEFRONTS 5   �  1      WAVEFRONTCUMSUM%SIZE+WAVEFRONTS=SIZE 3   �  0      WAVEFRONTCUMSUM%ABS+WAVEFRONTS=ABS 1      C   e   WAVEFRONTCUMSUM%BEAM1+WAVEFRONTS 1   b   C   e   WAVEFRONTCUMSUM%BEAM2+WAVEFRONTS *   �   �       CLEARWAVEFRONT+WAVEFRONTS @   G!  7      CLEARWAVEFRONT%ASSOCIATED+WAVEFRONTS=ASSOCIATED /   ~!  C   e   CLEARWAVEFRONT%BEAM+WAVEFRONTS $   �!  N      CLIPPOLY+WAVEFRONTS *   #  o      CLIPPOLY%POLYGON+POLYGONS B   ~#  u   %   CLIPPOLY%POLYGON%NUMVERTICES+POLYGONS=NUMVERTICES 8   �#  �   %   CLIPPOLY%POLYGON%VERTEX+POLYGONS=VERTEX ,   �$  0      CLIPPOLY%ABS+WAVEFRONTS=ABS 4   �$  4      CLIPPOLY%PRESENT+WAVEFRONTS=PRESENT 4   %  4      CLIPPOLY%EPSILON+WAVEFRONTS=EPSILON ,   ?%  0      CLIPPOLY%COS+WAVEFRONTS=COS ,   o%  0      CLIPPOLY%SIN+WAVEFRONTS=SIN )   �%  C   e   CLIPPOLY%BEAM+WAVEFRONTS )   �%  J   e   CLIPPOLY%POLY+WAVEFRONTS '   ,&  0   e   CLIPPOLY%XC+WAVEFRONTS '   \&  0   e   CLIPPOLY%YC+WAVEFRONTS *   �&  0   e   CLIPPOLY%ANGLE+WAVEFRONTS )   �&  0   e   CLIPPOLY%OBSC+WAVEFRONTS .   �&  0   e   CLIPPOLY%TILTANGLE+WAVEFRONTS 4   '  0   e   CLIPPOLY%TILTORIENTATION+WAVEFRONTS /   L'  �       WAVEFRONTWAVELENGTH+WAVEFRONTS 4   �'  C   e   WAVEFRONTWAVELENGTH%BEAM+WAVEFRONTS )   (  �       GENTOPHATBEAM+WAVEFRONTS 9   )  4      GENTOPHATBEAM%PRESENT+WAVEFRONTS=PRESENT 0   <)  C   e   GENTOPHATBEAM%TOPHAT+WAVEFRONTS ,   )  0   e   GENTOPHATBEAM%NX+WAVEFRONTS ,   �)  0   e   GENTOPHATBEAM%NY+WAVEFRONTS ,   �)  0   e   GENTOPHATBEAM%DX+WAVEFRONTS ,   *  0   e   GENTOPHATBEAM%DY+WAVEFRONTS 0   ?*  0   e   GENTOPHATBEAM%LAMBDA+WAVEFRONTS /   o*  0   e   GENTOPHATBEAM%DBEAM+WAVEFRONTS ,   �*  0   e   GENTOPHATBEAM%XC+WAVEFRONTS ,   �*  0   e   GENTOPHATBEAM%YC+WAVEFRONTS /   �*  0   e   GENTOPHATBEAM%POWER+WAVEFRONTS $   /+        CLIPCIRC+WAVEFRONTS ,   J,  0      CLIPCIRC%ABS+WAVEFRONTS=ABS 4   z,  4      CLIPCIRC%PRESENT+WAVEFRONTS=PRESENT 4   �,  4      CLIPCIRC%EPSILON+WAVEFRONTS=EPSILON ,   �,  0      CLIPCIRC%COS+WAVEFRONTS=COS )   -  C   e   CLIPCIRC%BEAM+WAVEFRONTS )   U-  0   e   CLIPCIRC%DIAM+WAVEFRONTS '   �-  0   e   CLIPCIRC%XC+WAVEFRONTS '   �-  0   e   CLIPCIRC%YC+WAVEFRONTS )   �-  0   e   CLIPCIRC%OBSC+WAVEFRONTS .   .  0   e   CLIPCIRC%TILTANGLE+WAVEFRONTS 4   E.  0   e   CLIPCIRC%TILTORIENTATION+WAVEFRONTS .   u.  �       APPLYARBITRARYMASK+WAVEFRONTS 8   M/  1      APPLYARBITRARYMASK%SIZE+WAVEFRONTS=SIZE >   ~/  4      APPLYARBITRARYMASK%PRESENT+WAVEFRONTS=PRESENT 3   �/  C   e   APPLYARBITRARYMASK%BEAM+WAVEFRONTS 3   �/  C   e   APPLYARBITRARYMASK%MASK+WAVEFRONTS 7   80  0   e   APPLYARBITRARYMASK%NEGATIVE+WAVEFRONTS    h0  �       FOCUSLENS    #1  C   a   FOCUSLENS%BEAM &   f1  0   a   FOCUSLENS%FOCUSLENS_F &   �1  0   a   FOCUSLENS%FOCUSLENS_D     �1  0   a   FOCUSLENS%DXOUT     �1  0   a   FOCUSLENS%DYOUT    &2        APPLYMASK "   =3  4      APPLYMASK%PRESENT    q3  C   a   APPLYMASK%BEAM     �3  0   a   APPLYMASK%LEN_R     �3  0   a   APPLYMASK%LEN_T !   4  0   a   APPLYMASK%OFFSET $   D4  8   a   APPLYMASK%DIRECTION #   |4  0   a   APPLYMASK%ROTANGLE "   �4  0   a   APPLYMASK%X_ALIGN "   �4  0   a   APPLYMASK%Y_ALIGN $   5  0   a   APPLYMASK%TILTANGLE *   <5  0   a   APPLYMASK%TILTORIENTATION     l5        APPLYMASKROTATE (   �6  4      APPLYMASKROTATE%PRESENT %   �6  C   a   APPLYMASKROTATE%BEAM &    7  0   a   APPLYMASKROTATE%LEN_R &   07  0   a   APPLYMASKROTATE%LEN_T '   `7  0   a   APPLYMASKROTATE%OFFSET )   �7  0   a   APPLYMASKROTATE%ROTANGLE *   �7  8   a   APPLYMASKROTATE%DIRECTION (   �7  0   a   APPLYMASKROTATE%X_ALIGN (   (8  0   a   APPLYMASKROTATE%Y_ALIGN *   X8  0   a   APPLYMASKROTATE%TILTANGLE 0   �8  0   a   APPLYMASKROTATE%TILTORIENTATION !   �8  �       APPLYMASKGENERAL &   g9  C   a   APPLYMASKGENERAL%BEAM '   �9  0   a   APPLYMASKGENERAL%LEN_X '   �9  0   a   APPLYMASKGENERAL%LEN_Y '   
:  0   a   APPLYMASKGENERAL%OFF_X '   ::  0   a   APPLYMASKGENERAL%OFF_Y    j:  �       APPLYMASKX     ;  1      APPLYMASKX%SQRT     K;  C   a   APPLYMASKX%BEAM "   �;  0   a   APPLYMASKX%RADIUS "   �;  0   a   APPLYMASKX%LENGTH "   �;  �       APPLYMASKXROUNDED '   �<  1      APPLYMASKXROUNDED%SQRT '   �<  C   a   APPLYMASKXROUNDED%BEAM )   )=  0   a   APPLYMASKXROUNDED%RADIUS )   Y=  0   a   APPLYMASKXROUNDED%LENGTH -   �=  0   a   APPLYMASKXROUNDED%CORNER_RAD    �=  �       APPLYMASKPOLY #   Z>  C   a   APPLYMASKPOLY%BEAM '   �>  x   a   APPLYMASKPOLY%VERTICES !   ?  0   a   APPLYMASKPOLY%XC !   E?  0   a   APPLYMASKPOLY%YC 